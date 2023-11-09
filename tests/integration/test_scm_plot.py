import os

from ciceroscm import CICEROSCM


def check_output(
    output_dir,
    expected_output_dir,
    files=[
        "conc_1.png",
        "conc_2.png",
        "em_1.png",
        "em_2.png",
        "em_3.png",
        "forc_1.png",
        "forc_2.png",
        "forc_3.png",
        "ohc.png",
        "rib.png",
        "temp.png",
    ],
):
    for filename in files:
        file_expected = os.path.join(expected_output_dir, "plots", filename)
        assert os.path.exists(file_expected.replace(expected_output_dir, output_dir))


def test_plot(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )

    outdir = str(tmpdir)

    cscm._run({"output_folder": outdir}, make_plot=True)

    check_output(outdir, os.path.join(test_data_dir, "ssp245_emis"))
