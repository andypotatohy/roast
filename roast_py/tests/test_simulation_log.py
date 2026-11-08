"""Tests for roast_py.simulation_log (writeRoastLog.m + isNewOptions.m) and
for roast()'s step checks: a repeated call reuses the tag and skips every
step, an interrupted one resumes, and redoing a step redoes the ones after
it. The heavy steps are replaced by fakes that count their calls."""

import json
import sys
import types
from datetime import datetime

import nibabel as nib
import numpy as np
import pytest

from roast_py import simulation_log as sl


def _options(**changes):
    opt = {
        "subj": "/data/subj.nii",
        "model_mri": "/data/subj.nii",
        "simulation_tag": None,
        "recipe": {"Fp1": 1.0, "P4": -1.0},
        "config_txt": "Fp1 (1 mA), P4 (-1 mA)",
        "cap_type": "1010",
        "electrodes": [
            {"elec_type": "disc", "elec_size": [6.0, 2.0], "elec_ori": None},
            {"elec_type": "pad", "elec_size": [50.0, 30.0, 3.0], "elec_ori": "lr"},
        ],
        "conductivities": {"white": 0.126, "gray": 0.276, "csf": 1.65, "bone": 0.01, "skin": 0.465,
                           "air": 2.5e-14, "gel": [0.3, 0.3], "electrode": [5.9e7, 5.9e7]},
        "mesh_options": {"radbound": 5.0, "angbound": 30.0, "distbound": 0.3, "reratio": 3.0, "maxvol": 10.0},
        "zero_padding": 0,
        "resampling": False,
        "segmentation": "multiaxial",
        "t2": None,
        "is_non_ras": False,
        "mri2mni": np.eye(4).tolist(),
        "Affine": np.eye(4).tolist(),
    }
    opt.update(changes)
    return opt


def test_timestamp_tag_is_matlabs_datestr_30():
    assert sl.timestamp_tag(datetime(2026, 10, 25, 14, 30, 12)) == "20261025T143012"


def test_same_options_are_not_new():
    assert not sl.is_new_options(_options(), json.loads(json.dumps(_options())))
    # 10-20, 10-10 and 10-05 are one cap to ROAST; a non-pad's orientation doesn't count
    assert not sl.is_new_options(_options(cap_type="1020"), _options())
    electrodes = [dict(_options()["electrodes"][0], elec_ori="ap"), _options()["electrodes"][1]]
    assert not sl.is_new_options(_options(electrodes=electrodes), _options())


@pytest.mark.parametrize(
    "changes",
    [
        {"recipe": {"Fp1": 0.5, "P4": -0.5}},
        {"recipe": {"P4": -1.0, "Fp1": 1.0}},
        {"cap_type": "biosemi"},
        {"t2": "/data/subj_T2_aligned.nii"},
        {"segmentation": "spm"},
        {"Affine": (np.eye(4) * 1.01).tolist()},
        {"mesh_options": {"radbound": 4.0, "angbound": 30.0, "distbound": 0.3, "reratio": 3.0, "maxvol": 10.0}},
        {"resampling": True},
        {"zero_padding": 10},
    ],
)
def test_changed_options_are_new(changes):
    assert sl.is_new_options(_options(**changes), _options())


def test_electrode_and_conductivity_changes_are_new():
    base = _options()
    sizes = [dict(base["electrodes"][0], elec_size=[8.0, 2.0]), base["electrodes"][1]]
    assert sl.is_new_options(_options(electrodes=sizes), base)
    oris = [base["electrodes"][0], dict(base["electrodes"][1], elec_ori=[0.71, 0.71, 0.0])]
    assert sl.is_new_options(_options(electrodes=oris), base)  # keyword vs vector, as in MATLAB
    cond = dict(base["conductivities"], gel=[1.0, 0.3])
    assert sl.is_new_options(_options(conductivities=cond), base)
    # options saved by an older roast_py, missing fields
    assert sl.is_new_options(base, {k: v for k, v in base.items() if k != "segmentation"})


def test_register_simulation_tags_reuses_and_refuses(tmp_path):
    work = str(tmp_path)
    first, reused = sl.register_simulation(_options(), work, "subj")
    assert not reused and len(first["simulation_tag"]) == 15  # yyyymmddTHHMMSS
    saved = tmp_path / f"subj_{first['simulation_tag']}_roastOptions.json"
    assert json.loads(saved.read_text())["simulation_tag"] == first["simulation_tag"]
    log = (tmp_path / "subj_roastLog").read_text()
    assert log.startswith(f"{first['simulation_tag']}:\nrecipe:\tFp1 (1 mA), P4 (-1 mA)\n")
    assert "elecOri:\tN/A; lr\n" in log and "multiaxial:\ton\n" in log

    again, reused = sl.register_simulation(_options(), work, "subj")
    assert reused and again["simulation_tag"] == first["simulation_tag"]
    with pytest.warns(UserWarning, match="will be ignored"):
        named, reused = sl.register_simulation(_options(simulation_tag="mine"), work, "subj")
    assert reused and named["simulation_tag"] == first["simulation_tag"]

    other, reused = sl.register_simulation(_options(simulation_tag="mine", zero_padding=10), work, "subj")
    assert not reused and other["simulation_tag"] == "mine"
    with pytest.raises(ValueError, match="Please use a new tag"):
        sl.register_simulation(_options(simulation_tag="mine", zero_padding=20), work, "subj")
    assert log in (tmp_path / "subj_roastLog").read_text()  # appended, not rewritten

    # another MRI whose name starts with "subj_" isn't mistaken for this one
    sl.register_simulation(_options(subj="/data/subj_1mm.nii", simulation_tag="x"), work, "subj_1mm")
    assert [o["simulation_tag"] for o in sl.saved_simulations(work, "subj", "/data/subj.nii")] == [
        first["simulation_tag"], "mine",
    ]
    assert sl.latest_simulation(work, "subj", "/data/subj.nii")["simulation_tag"] == "mine"


def test_two_untagged_runs_in_one_second_get_distinct_tags(tmp_path, monkeypatch):
    monkeypatch.setattr(sl, "timestamp_tag", lambda now=None: "20261025T143012")
    a, _ = sl.register_simulation(_options(), str(tmp_path), "subj")
    b, _ = sl.register_simulation(_options(zero_padding=5), str(tmp_path), "subj")
    assert (a["simulation_tag"], b["simulation_tag"]) == ("20261025T143012", "20261025T143012-2")


@pytest.mark.parametrize("tag", ["", "a/b", "x:y"])
def test_bad_tags_are_refused(tag):
    with pytest.raises(ValueError, match="Invalid simulation tag"):
        sl.check_tag(tag)


# -- roast()'s step checks ----------------------------------------------------


class _FakePipeline:
    """Stand-ins for the six steps, each counting its calls and writing the
    files the real one writes."""

    def __init__(self, monkeypatch, shape=(12, 12, 12)):
        import roast_py.fem.solve as fem_solve
        import roast_py.geometry.placement as placement
        import roast_py.meshing.cgal_mesher as mesher
        import roast_py.registration as registration
        from roast_py.registration import Registration, registration_path

        self.calls = {k: 0 for k in ("segment", "register", "place", "mesh", "solve", "post")}
        self.shape = shape

        def segment(t1, out_path=None, model_dir=None):
            self.calls["segment"] += 1
            labels = np.zeros(shape, np.uint8)
            labels[2:-2, 2:-2, 2:-2] = 5
            nib.save(nib.Nifti1Image(labels, nib.load(t1).affine), out_path)
            return out_path

        def run_niftyreg(mri, out_dir=None, bin_path=None):
            self.calls["register"] += 1
            reg = Registration(np.eye(4), nib.load(mri).affine, np.eye(4))
            reg.save(registration_path(mri, out_dir))
            return reg

        def electrode_placement(labels, *args, **kwargs):
            self.calls["place"] += 1
            elec = np.zeros(shape, np.uint8)
            elec[0, 0, :2] = [1, 2]
            return elec, elec.copy()

        def mesh_by_iso2mesh(labels, elec, gel, voxel_size, work_dir, out_path, bin_path=None, **opts):
            self.calls["mesh"] += 1
            with open(out_path, "w") as f:
                f.write("$MeshFormat\n")
            return np.zeros((4, 4)), np.array([[1, 2, 3, 4, 1]]), np.zeros((1, 4))

        def solve(work_dir, stem, node, elem, names, current, sigma, bin_path=None):
            self.calls["solve"] += 1
            for path, n in zip(fem_solve.solution_paths(work_dir, stem), (1, 3)):
                with open(path, "w") as f:
                    f.write("2\n" + "".join(f"{i} " + " ".join(["0.5"] * n) + "\n" for i in (1, 2)))

        def postprocess(work_dir, stem, node, voxel_size, grid_shape):
            self.calls["post"] += 1
            v = np.full(grid_shape, 1.0)
            return v, np.zeros((*grid_shape, 3)), v * 2

        fake_multiaxial = types.ModuleType("roast_py.segmentation.multiaxial")
        fake_multiaxial.segment = segment  # no TensorFlow import
        monkeypatch.setitem(sys.modules, "roast_py.segmentation.multiaxial", fake_multiaxial)
        monkeypatch.setattr(registration, "run_niftyreg", run_niftyreg)
        monkeypatch.setattr(placement, "electrode_placement", electrode_placement)
        monkeypatch.setattr(mesher, "mesh_by_iso2mesh", mesh_by_iso2mesh)
        monkeypatch.setattr(fem_solve, "solve", solve)
        monkeypatch.setattr(fem_solve, "postprocess", postprocess)

    def counts(self):
        out = dict(self.calls)
        for k in self.calls:
            self.calls[k] = 0
        return out


@pytest.fixture
def fake_run(tmp_path, monkeypatch):
    fake = _FakePipeline(monkeypatch)
    subj = tmp_path / "head.nii"
    nib.save(nib.Nifti1Image(np.ones(fake.shape, np.float32), np.eye(4)), subj)

    def run(**kwargs):
        from roast_py import roast

        return roast(str(subj), install_missing=False, visualize=False, **kwargs)

    return fake, run, tmp_path


ALL = {"segment": 1, "register": 1, "place": 1, "mesh": 1, "solve": 1, "post": 1}
NONE = dict.fromkeys(ALL, 0)


def test_roast_skips_what_exists_and_reuses_the_tag(fake_run):
    fake, run, tmp = fake_run
    first = run()
    assert fake.counts() == ALL and len(first.sim_tag) == 15
    tag = first.sim_tag
    for suffix in ("_roastOptions.json", "_mask_elec.nii", "_mask_gel.nii", "_mesh.npz", ".msh", "_v.pos",
                   "_e.pos", "_v.nii", "_e.nii", "_emag.nii"):
        assert (tmp / f"head_{tag}{suffix}").exists(), suffix
    assert not list(tmp.glob(".partial_*"))

    again = run()
    assert fake.counts() == NONE and again.sim_tag == tag
    assert np.array_equal(again.ef_mag, first.ef_mag) and again.ef_mag.dtype == first.ef_mag.dtype

    # a different recipe is a new simulation of the same head: steps 1-2 are shared
    other = run(recipe={"Fp1": 1, "Cz": -1}, simulation_tag="cz")
    assert fake.counts() == dict(ALL, segment=0, register=0) and other.sim_tag == "cz"
    assert len(list(tmp.glob("head_*_roastOptions.json"))) == 2
    assert (tmp / "head_roastLog").read_text().count("recipe:\t") == 2


def test_roast_resumes_and_redoes_the_steps_after_a_redone_one(fake_run):
    fake, run, tmp = fake_run
    tag = run(simulation_tag="t").sim_tag
    fake.counts()

    (tmp / f"head_{tag}_v.nii").unlink()  # interrupted while saving the results
    run(simulation_tag="t")
    assert fake.counts() == dict(NONE, post=1)

    with open(tmp / f"head_{tag}_e.pos", "w") as f:  # getDP cut off mid-file
        f.write("2\n1 0.5 0.5 0.5\n")
    run(simulation_tag="t")
    assert fake.counts() == dict(NONE, solve=1, post=1)

    (tmp / f"head_{tag}_mask_elec.nii").unlink()  # redo the placement: and so the mesh, solve, results
    run(simulation_tag="t")
    assert fake.counts() == dict(NONE, place=1, mesh=1, solve=1, post=1)


def test_roast_refuses_a_used_tag_for_new_options(fake_run):
    fake, run, tmp = fake_run
    run(simulation_tag="t")
    fake.counts()
    with pytest.raises(ValueError, match="Please use a new tag"):
        run(simulation_tag="t", zero_padding=2)
    assert fake.counts()["place"] == 0
    with pytest.raises(ValueError, match="Invalid simulation tag"):
        run(simulation_tag="a/b")
