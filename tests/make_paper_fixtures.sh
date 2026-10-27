#!/bin/bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Utkarsh Singh
# Fixtures from the paper's benchmark runs (VASP 6.5.1), which are not part of this repository:
#     PAPER_RUNS=/path/to/the/runs [PYTHON=python3] ./make_paper_fixtures.sh
# e1 = stock IBRION=6 ("Standard" column of Table IV), e2 = LMINIMAL_ELASTIC v1 ("This work").
set -e
PY=${PYTHON:-python}
R=${PAPER_RUNS:?set PAPER_RUNS to the directory of the paper's benchmark runs}
H=$(cd "$(dirname "$0")" && pwd)
O=$H/data_paper
# run from inside $R so that the fixtures record the run directory relative to it
mk() { (cd "$R" && $PY "$H/make_vasp_fixture.py" "$1" "$O/$2.json" --label "$3" ${4:+--paper "$4"} ${5:+--paper-upper "$5"}); }
mk Al/e1 al_e1 "fcc Al, Table IV Standard" "c11=114.2,c12=59.8,c44=34.2"
mk Al/e2 al_e2 "fcc Al, Table IV This work" "c11=115.3,c12=59.1,c44=34.3"
mk Mg/e1/s2 mg_e1 "hcp Mg, Table IV Standard" "c11=59.8,c12=26.7,c13=20.6,c33=62.5,c44=16.9"
mk Mg/e2/s2 mg_e2 "hcp Mg, Table IV This work" "c11=60.2,c12=26.6,c13=20.5,c33=62.5,c44=17.0"
mk SiO2/e1 sio2_e1 "alpha-SiO2 P3_221, Table IV Standard" "c11=79.2,c12=-3.0,c13=6.9,c14=18.3,c33=87.0,c44=51.2"
mk SiO2/e2 sio2_e2 "alpha-SiO2 P3_221, Table IV This work" "c11=79.7,c12=-2.7,c13=7.4,c14=18.1,c33=87.7,c44=51.6"
mk "CaMg(Co3)2/e1" dolomite_e1 "dolomite R-3, Table IV Standard" "c11=191.3,c12=65.9,c13=54.6,c14=16.7,c15=10.9,c33=106.5,c44=37.7"
mk "CaMg(Co3)2/e2" dolomite_e2 "dolomite R-3, Table IV This work" "c11=192.6,c12=67.9,c13=55.4,c14=16.7,c15=10.9,c33=106.3,c44=38.1"
mk TiO2/e1 tio2_e1 "rutile TiO2, Table IV Standard" "c11=256.6,c12=170.2,c13=147.6,c33=467.9,c44=113.9,c66=211.6"
mk TiO2/e2 tio2_e2 "rutile TiO2, Table IV This work" "c11=256.8,c12=170.7,c13=147.8,c33=467.9,c44=113.9,c66=211.9"
mk CaMoO4/e1 camoo4_e1 "CaMoO4 I4_1/a, Table IV Standard" "c11=131.8,c12=52.6,c13=46.4,c16=-9.9,c33=114.3,c44=31.7,c66=38.7"
mk CaMoO4/e2 camoo4_e2 "CaMoO4 I4_1/a, Table IV This work" "c11=131.9,c12=52.7,c13=46.4,c16=-9.9,c33=114.2,c44=31.8,c66=38.7"
# TiSi2 and NaAlSi3O8: the table is in another axis frame -> full upper triangle, compared by eigenvalues
mk TiSi2/e1 tisi2_e1 "TiSi2 Fddd, Table IV Standard (paper frame: x and y exchanged)" "" "312.2,27.2,26.7,0,0,0,313.3,82.3,0,0,0,399.6,0,0,0,106.6,0,0,73.1,0,115.1"
mk TiSi2/e2 tisi2_e2 "TiSi2 Fddd, Table IV This work (paper frame: x and y exchanged)" "" "312.2,27.3,26.7,0,0,0,313.4,82.5,0,0,0,399.6,0,0,0,106.3,0,0,72.9,0,115.3"
mk NaAlSi3O8/e1 albite_e1 "NaAlSi3O8 P-1, Table IV Standard (paper in a rotated frame)" "" "63.1,27.6,28.3,5.7,-5.1,-0.1,162.1,11.8,-3.6,-7.8,-4.0,163.2,-5.5,5.0,-8.9,24.2,-1.4,-6.7,26.9,1.1,33.6"
mk NaAlSi3O8/e2 albite_e2 "NaAlSi3O8 P-1, Table IV This work (paper in a rotated frame)" "" "63.1,27.6,28.3,5.7,-5.2,-0.1,162.2,11.8,-3.6,-7.8,-4.0,163.2,-5.5,5.0,-8.9,24.2,-1.4,-6.7,27.0,1.1,33.6"
# monoclinic: sibling run of the Table IV beta-Ga2O3 (slightly different settings; eigenvalues within ~2%)
mk beta-Ga2O3/e1/t1 ga2o3_e1 "beta-Ga2O3 C2/m, stock (sibling of the Table IV run)"
mk beta-Ga2O3/e2/t1 ga2o3_e2 "beta-Ga2O3 C2/m, LMINIMAL_ELASTIC v1 (sibling of the Table IV run)"
# extra materials of the same campaign
mk CaWO4/e1 cawo4_e1 "CaWO4 I4_1/a, stock"
mk CaWO4/e2 cawo4_e2 "CaWO4 I4_1/a, LMINIMAL_ELASTIC v1"
mk Bi/e1 bi_e1 "Bi R-3m, stock"
mk Bi/e2 bi_e2 "Bi R-3m, LMINIMAL_ELASTIC v1"
mk BaTiO3/e1 batio3_e1 "BaTiO3 P4mm, stock"
mk BaTiO3/e2 batio3_e2 "BaTiO3 P4mm, LMINIMAL_ELASTIC v1"
