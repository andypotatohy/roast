"""Ports the .pro (getDP "problem" script) generation from solveByGetDP.m:
defines the electrostatics/Laplace problem (Region groups, conductivities,
Neumann current-density boundary conditions, the weak-form Galerkin
equation, and what to write out) that getDP actually solves.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field


@dataclass
class Conductivities:
    """mirrors ROAST's `conductivities` option struct. `gel`/`electrode`
    are per-electrode (indexed the same way as `current`/`indUse` below),
    matching ROAST's support for per-electrode conductivity overrides."""

    white: float = 0.126
    gray: float = 0.276
    csf: float = 1.65
    bone: float = 0.01
    skin: float = 0.465
    air: float = 2.5e-14
    gel: list[float] = field(default_factory=list)
    electrode: list[float] = field(default_factory=list)


def write_pro_file(
    pro_path: str,
    current: list[float],
    sigma: Conductivities,
    area_elec_needed,
    ind_use: list[int],
    lf_tag: str = "",
) -> None:
    """Ports solveByGetDP.m's .pro generation.

    `current` and `sigma.gel`/`sigma.electrode` are indexed like ROAST's
    MATLAB arrays: `current[e]`/`sigma.gel[e]` for the (1-based, via
    `ind_use`) electrode `e`. `ind_use` (1-based) selects which
    electrodes' boundary conditions actually get written into this
    Resolution -- for a full roast() solve this is every electrode; for
    lead-field generation (`lf_tag` nonempty, one call per electrode)
    it's a single electrode at a time, matching solveByGetDP.m's own
    `indUse`/`LFtag` parameters.
    """
    num_tissue = 6
    num_elec = len(area_elec_needed)

    lines: list[str] = []

    def w(s: str = "") -> None:
        lines.append(s)

    w("Group {")
    w()
    w("white = Region[1];")
    w("gray = Region[2];")
    w("csf = Region[3];")
    w("bone = Region[4];")
    w("skin = Region[5];")
    w("air = Region[6];")
    for i in ind_use:
        w(f"gel{i} = Region[{num_tissue + i}];")
    for i in ind_use:
        w(f"elec{i} = Region[{num_tissue + num_elec + i}];")

    gel_str = "".join(f"gel{i}, " for i in ind_use)
    elec_str = "".join(f"elec{i}, " for i in ind_use)
    used_elec_str = ""
    for i in ind_use:
        used_elec_str += f"usedElec{i}, "
        w(f"usedElec{i} = Region[{num_tissue + 2 * num_elec + i}];")

    w(f"DomainC = Region[{{white, gray, csf, bone, skin, air, {gel_str}{elec_str[:-2]}}}];")
    w()
    w(f"AllDomain = Region[{{white, gray, csf, bone, skin, air, {gel_str}{elec_str}{used_elec_str[:-2]}}}];")
    w()
    w("}")
    w()

    w("Function {")
    w()
    w(f"sigma[white] = {sigma.white:g};")
    w(f"sigma[gray] = {sigma.gray:g};")
    w(f"sigma[csf] = {sigma.csf:g};")
    w(f"sigma[bone] = {sigma.bone:g};")
    w(f"sigma[skin] = {sigma.skin:g};")
    w(f"sigma[air] = {sigma.air:g};")
    for i in ind_use:
        w(f"sigma[gel{i}] = {sigma.gel[i - 1]:g};")
    for i in ind_use:
        w(f"sigma[elec{i}] = {sigma.electrode[i - 1]:g};")

    for i in ind_use:
        du_dn = 1000 * current[i - 1] / area_elec_needed[i - 1]
        w(f"du_dn{i}[] = {du_dn:g};")

    w()
    w("}")
    w()

    w("Jacobian {")
    w("  { Name Vol ;")
    w("    Case {")
    w("      { Region All ; Jacobian Vol ; }")
    w("    }")
    w("  }")
    w("  { Name Sur ;")
    w("    Case {")
    w("      { Region All ; Jacobian Sur ; }")
    w("    }")
    w("  }")
    w("}")
    w()

    w("Integration {")
    w("  { Name GradGrad ;")
    w("    Case { {Type Gauss ;")
    w("            Case { { GeoElement Triangle    ; NumberOfPoints  3 ; }")
    w("                   { GeoElement Quadrangle  ; NumberOfPoints  4 ; }")
    w("                   { GeoElement Tetrahedron ; NumberOfPoints  4 ; }")
    w("                   { GeoElement Hexahedron  ; NumberOfPoints  6 ; }")
    w("                   { GeoElement Prism       ; NumberOfPoints  9 ; } }")
    w("           }")
    w("         }")
    w("  }")
    w("}")
    w()

    w("FunctionSpace {")
    w("  { Name Hgrad_v_Ele; Type Form0;")
    w("    BasisFunction {")
    w("      { Name sn; NameOfCoef vn; Function BF_Node;")
    w("        Support AllDomain; Entity NodesOf[ All ]; }")
    w("    }")
    w("  }")
    w("}")
    w()

    w("Formulation {")
    w("  { Name Electrostatics_v; Type FemEquation;")
    w("    Quantity {")
    w("      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }")
    w("    }")
    w("    Equation {")
    w("      Galerkin { [ sigma[] * Dof{d v} , {d v} ]; In DomainC; ")
    w("                 Jacobian Vol; Integration GradGrad; }")
    w()
    for i in ind_use:
        w(f"      Galerkin{{ [ -du_dn{i}[], {{v}} ]; In usedElec{i};")
        w("                 Jacobian Sur; Integration GradGrad;}")
    w("    }")
    w("  }")
    w("}")
    w()

    w("Resolution {")
    w("  { Name EleSta_v;")
    w("    System {")
    w("      { Name Sys_Ele; NameOfFormulation Electrostatics_v; }")
    w("    }")
    w("    Operation { ")
    w("      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];")
    w("    }")
    w("  }")
    w("}")
    w()

    w("PostProcessing {")
    w("  { Name EleSta_v; NameOfFormulation Electrostatics_v;")
    w("    Quantity {")
    w("      { Name v; ")
    w("        Value { ")
    w("          Local { [ {v} ]; In AllDomain; Jacobian Vol; } ")
    w("        }")
    w("      }")
    w("      { Name e; ")
    w("        Value { ")
    w("          Local { [ -{d v} ]; In AllDomain; Jacobian Vol; }")
    w("        }")
    w("      }")
    w("    }")
    w("  }")
    w("}")
    w()

    # getDP's Print directives below use a bare filename (matching
    # solveByGetDP.m, which does the same) -- the caller must run getdp
    # with cwd set to this .pro file's directory for the output .pos files
    # to land next to it (see getdp_runner.py).
    subj_name = os.path.splitext(os.path.basename(pro_path))[0]

    w("PostOperation {")
    w()
    w("{ Name Map; NameOfPostProcessing EleSta_v;")
    w("   Operation {")
    if not lf_tag:
        w(f'     Print [ v, OnElementsOf DomainC, File "{subj_name}_v.pos", Format NodeTable ];')
    w(f'     Print [ e, OnElementsOf DomainC, Smoothing, File "{subj_name}_e{lf_tag}.pos", Format NodeTable ];')
    w("   }")
    w("}")
    w()
    w("}")

    with open(pro_path, "w") as f:
        f.write("\n".join(lines) + "\n")
