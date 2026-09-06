"""Only the batch-job text generation is testable without a real SPM12
standalone + MATLAB Runtime install (not available in this environment --
see spm_standalone.py's module docstring)."""

from roast_py.segmentation.spm_standalone import _TISSUE_NGAUS, write_batch_job


def test_write_batch_job_matches_start_seg_structure(tmp_path):
    t1 = tmp_path / "subject1.nii"
    t1.write_bytes(b"")  # content doesn't matter, only the path is embedded
    template = tmp_path / "eTPM.nii"
    template.write_bytes(b"")

    job_path = tmp_path / "job.m"
    write_batch_job(str(t1), str(job_path), template_path=str(template))

    script = job_path.read_text()

    assert f"{{'{t1},1'}}" in script
    assert "matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;" in script
    assert "matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];" in script

    # 6 tissue classes, each referencing the right TPM volume index and the
    # ngaus values start_seg.m uses (tied to the bundled eTPM.nii atlas).
    for i, ngaus in enumerate(_TISSUE_NGAUS, start=1):
        assert f"tissue({i}).tpm = {{'{template},{i}'}};" in script
        assert f"tissue({i}).ngaus = {ngaus};" in script
        assert f"tissue({i}).native = [1 0];" in script

    assert "spm_jobman('initcfg');" in script
    assert "spm_jobman('run', matlabbatch);" in script


def test_write_batch_job_native_plus_normalized(tmp_path):
    t1 = tmp_path / "subject1.nii"
    t1.write_bytes(b"")
    template = tmp_path / "eTPM.nii"
    template.write_bytes(b"")
    job_path = tmp_path / "job.m"

    write_batch_job(str(t1), str(job_path), template_path=str(template), native=False)
    script = job_path.read_text()
    assert "tissue(1).native = [1 1];" in script
