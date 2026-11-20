from pathlib import Path

from openalea.simbal.tools import loading_parameters_from_file


def test_load_parameters_from_file():

    BASE_DIR = Path(__file__).resolve().parent
    parameters_filepath = BASE_DIR / "data" / "simbal_parameters.xlsx"
    outputs_path = BASE_DIR / "data" / "outputs"

    parameters = loading_parameters_from_file(parameters_filepath = parameters_filepath, outputs_path=outputs_path)

    assert parameters['simulation_period_in_days'] == 60