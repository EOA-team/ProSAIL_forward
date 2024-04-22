# ProSAIL_forward

This repository allows to run the ProSAIL RTM in forward mode and generate Look-Up Tables (LUTs).
It provides default data and parametrization for running the model for winter wheat in Switzerland and for the Sentinel-2 sensor.


## Installation

```
git clone https://github.com/EOA-team/ProSAIL_forward.git
cd ProSAIL_forward
python -m venv my_venv
source my_venv/bin/activate
pip install -r requirements.txt
```

## Running the code

To run the RTM, everything is set up in the `simulate_S2_spectra_soil.py` script. 

First define some paths and parameters:

| Parameter        | Description                                                                                                   |
|------------------|---------------------------------------------------------------------------------------------------------------|
| `out_dir`        | The output directory where results are saved.                                                                 |
| `rtm_lut_config` | Dictionary containing number of simulations and sensor (Sentinel-2A and 2B must be run separately!).          |
| `lut_params_dir` | Directory where the RTM variables are saved. The code will look for 'prosail_danner-etal_switzerland.csv' by default. If named differently, you can edit the code to point to the correct file. |
| `soil_path`      | Optional. If you want the RTM to use custom soil samples, then pass the path to a .pkl dataframe containing each spectra in a row (2101 columns with reflectance values between 400 and 2500nm). If none provided, it will use default soil spectra. |
| `traits`         | Traits to include among LAI (`lai`), chlorophyll (`cab`), carotenoids (`car`), canopy chlorophyll content (`ccc`). Default is all four. |



You can run it with the following command:

```
python simulate_S2_spectra_soil.py
```

The result will be a pickled dataframe containing a simulation per row. The columns will be the leaf/canopy parameters and the Sentinel2 bands' reflectance values:
```
import pandas as pd
file_path = 'path_to/your_results.pkl'
df = pd.read_pickle(file_path)
print(df.head()) 
```

## Credits

The code provided here is built on existing code
- `rtm_inv`: A Python-backend for radiative transfer model inversion for crop trait retrieval. [Code](https://github.com/EOA-team/rtm_inv)
- `sentinel2_crop_traits`: Sentinel-2 Crop Trait Retrieval Using Physiological and Phenological Priors from Field Phenotyping (Graf et al., 2023, RSE). [Code](https://github.com/EOA-team/sentinel2_crop_traits)



