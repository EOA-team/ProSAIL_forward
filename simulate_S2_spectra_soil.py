"""
Generate PROSAIL RTM simulations with soil spectra and a corresponding look up table.

@author Selene Ledain
"""

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional
from datetime import datetime
from rtm_inv.core.lookup_table import LookupTable, generate_lut, simulate_from_lut
import pickle
import numpy as np
import pandas as pd
import yaml


def get_logger():
  """
  Returns a logger object with stream and file handler
  """
  
  CURRENT_TIME: str = datetime.now().strftime("%Y%m%d-%H%M%S")
  LOGGER_NAME: str = "CropCovEO"
  LOG_FORMAT: str = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
  LOG_DIR: str = str(Path.home())  # ..versionadd:: 0.2.1
  LOG_FILE: str = os.path.join(LOG_DIR, f"{CURRENT_TIME}_{LOGGER_NAME}.log")
  LOGGING_LEVEL: int = logging.INFO

  # create file handler which logs even debug messages
  logger = logging.getLogger(LOGGER_NAME)
  logger.setLevel(LOGGING_LEVEL)
  
  fh: logging.FileHandler = logging.FileHandler(LOG_FILE)
  fh.setLevel(LOGGING_LEVEL)
  # create console handler with a higher log level
  ch: logging.StreamHandler = logging.StreamHandler()
  ch.setLevel(LOGGING_LEVEL)
  # create formatter and add it to the handlers
  formatter: logging.Formatter = logging.Formatter(LOG_FORMAT)
  fh.setFormatter(formatter)
  ch.setFormatter(formatter)
  # add the handlers to the logger
  logger.addHandler(fh)
  logger.addHandler(ch)

  return logger


def codistribute_vars(lut, lut_config):
    """ 
    Linearly codistribute the variables with LAI based of codsitrbution functions from a table

    :param lut: LookupTable object
    :param lut_config: dict
    """
    if lut_config['codistribution'] is not None:
      codist = pd.read_csv(lut_config['codistribution'])
      for i, row in codist.iterrows():
          param = row['param']
          Vmin0 = row['Vmin0']
          Vmax0 = row['Vmax0']
          VminLAI = row['Vmin(LAImax)']
          VmaxLAI = row['Vmax(LAImax)']
          lut._samples[param] = lut._samples[param].apply(lambda x: (x-Vmin0)/(Vmax0-Vmin0)*(VmaxLAI-VminLAI)+VminLAI)
 
    return lut


def generate_spectra(
    output_dir: Path,
    lut_params: Path,
    lut_config: Dict[str, Any],
    rtm_config: Dict[str, Any],
    traits: List[str]
  ) -> None:

  logger = get_logger()
  
  # run PROSAIL forward runs for the different parametrizations available
  logger.info('Starting PROSAIL runs')
  
  sensor_suffix = ''
  if rtm_config['sensor'] == 'Sentinel2A':
    sensor_suffix = '_S2A'  
  elif rtm_config['sensor'] == 'Sentinel2B':
    sensor_suffix = '_S2B'
  elif rtm_config['sensor'] == 'PlanetSuperDove':
    sensor_suffix = '_PL'
    
  pheno_phases = \
      lut_params.name.split('.csv')[0] + sensor_suffix

  # generate lookup-table
  trait_str = '-'.join(traits)
  fpath_lut = output_dir.joinpath(
    f'{pheno_phases}_{trait_str}_lut.pkl') 
  print('Output path', fpath_lut)

  # if LUT exists, continue, else generate it
  if not fpath_lut.exists():
    lut_inp = lut_config.copy()
    del lut_inp['codistribution']
    lut_inp['lut_params'] = lut_params
    lut = generate_lut(**lut_inp)
    if lut_config['codistribution'] is not None:
      lut = codistribute_vars(lut, lut_config) # codistribution
    lut = simulate_from_lut(lut, **rtm_config)

    # special case CCC (Canopy Chlorophyll Content) ->
    # this is not a direct RTM output
    if 'ccc' in traits:
      lut['ccc'] = lut['lai'] * lut['cab']
      # convert to g m-2 as this is the more common unit
      # ug -> g: factor 1e-6; cm2 -> m2: factor 1e-4
      lut['ccc'] *= 1e-2


    lut.dropna(inplace=True)

    # save LUT to file
    with open(fpath_lut, 'wb') as f:
      pickle.dump(lut, f, protocol=3)

  else:
    pass

  logger.info('Finished PROSAIL runs')


def generate_spectra_soil(
  output_dir: Path,
  lut_params: Path,
  lut_config: Dict[str, Any],
  rtm_config: Dict[str, Any],
  traits: List[str],
  soil_df: pd.DataFrame
  ) -> None:

  logger = get_logger()
  
  # run PROSAIL forward runs for the different parametrizations available
  logger.info('Starting PROSAIL runs')
  logger.info(rtm_config['sensor'])
  
  sensor_suffix = ''
  if rtm_config['sensor'] == 'Sentinel2A':
    sensor_suffix = '_S2A'  
  elif rtm_config['sensor'] == 'Sentinel2B':
    sensor_suffix = '_S2B'
  elif rtm_config['sensor'] == 'PlanetSuperDove':
    sensor_suffix = '_PL'
  elif rtm_config['sensor'] == 'Hyspex':
    sensor_suffix = '_hyspex'
  

  pheno_phases = \
      lut_params.name.split('.csv')[0] + sensor_suffix

  # generate lookup-table
  trait_str = '-'.join(traits)
  fpath_lut = output_dir.joinpath(
    f'{pheno_phases}_{trait_str}_lut.pkl') 
  print('Output path', fpath_lut)

  # if LUT exists, continue, else generate it
  if not fpath_lut.exists():
    # Generate LUTd
    lut_inp = lut_config.copy()
    del lut_inp['codistribution']
    lut_inp['lut_params'] = lut_params
    lut = generate_lut(**lut_inp)
    if lut_config['codistribution'] is not None:
      lut = codistribute_vars(lut, lut_config) # codistribution

    # Simulate with RTM
    rtm_inp = rtm_config.copy()
    # Create LUT subgroups of size lut_size/len(soil_df). 
    # Pass each subgroup with a soil spectra and simulate.
    # Concatenate all simulations
    sub_luts = get_random_subgroups(lut._samples, len(soil_df))
    lut_allsoils = []
    for i, (original_idx, sub_lut) in enumerate(sub_luts):
      logger.info(f'Simulating with soil spectra {i+1} of {len(soil_df)}')
      rtm_inp['rsoil0'] = None
      rtm_inp['soil_spectrum1'] = soil_df.iloc[i].values 
      rtm_inp['soil_spectrum2'] = np.zeros_like(soil_df.iloc[i].values) 
      sub_lut = dataframe_to_lookup_table(sub_lut, lut)
      lut_soilspectra = simulate_from_lut(sub_lut, **rtm_inp)
      lut_soilspectra.index = original_idx # keep original order
      lut_allsoils.append(lut_soilspectra)
    
    lut = pd.concat(lut_allsoils).sort_index()


    # special case CCC (Canopy Chlorophyll Content) ->
    # this is not a direct RTM output
    if 'ccc' in traits:
      lut['ccc'] = lut['lai'] * lut['cab']
      # convert to g m-2 as this is the more common unit
      # ug -> g: factor 1e-6; cm2 -> m2: factor 1e-4
      lut['ccc'] *= 1e-2

    # prepare LUT for model training
    # lut = lut[band_selection + traits].copy()
    lut.dropna(inplace=True)

 
    # save LUT to file
    with open(fpath_lut, 'wb+') as f:
      pickle.dump(lut, f)

  else:
    pass

  logger.info('Finished PROSAIL runs')


def get_random_subgroups(df, n):
    # Step 1: Shuffle the indices of the DataFrame
    indices = np.arange(len(df))
    np.random.shuffle(indices)
    
    # Step 2: Split the shuffled indices into `n` subgroups of roughly equal size
    subgroups = np.array_split(indices, n)
    
    # Step 3: Create subgroups of the DataFrame using these indices
    subgroups_dfs = [(indices, df.iloc[indices]) for indices in subgroups]
    
    return subgroups_dfs


def dataframe_to_lookup_table(df: pd.DataFrame, original_lut: LookupTable) -> LookupTable:
    """
    Converts a DataFrame to a LookupTable object.

    :param df: DataFrame to convert
    :param original_lut: Original LookupTable object to copy parameters from
    :return: LookupTable object
    """
    lut = LookupTable(original_lut._params_df)
    lut.samples = df.reset_index(drop=True).copy()
    return lut


def load_config(config_path: str) -> Dict:
  ''' 
  Load configuration file

  :param config_path: path to yaml file
  :returns: dictionary of parameters
  '''
  with open(config_path, "r") as config_file:
      config = yaml.safe_load(config_file)
  return config


if __name__ == '__main__':


  cwd = Path(__file__).parent.absolute()
  import os
  os.chdir(cwd)

  # Load configuration file 
  config_path = 'config.yaml'
  config = load_config(config_path)


  ########################
  # EXTRACT PARAMETERS AND PATHS FOR RUNNING

  # Store results/simulations
  out_dir = Path(config['out_dir'])
  out_dir.mkdir(exist_ok=True)

  # Prepare other params
  lut_params = Path(config['LUT'].pop('lut_params')) # TO DO
  fpath_srf = Path(config['RTM']['fpath_srf'])
  lut_config = config['LUT']
  rtm_config = config['RTM']
  soil_path = Path(config['soil_path']) if config['soil_path'] is not None else None 
  traits = config['traits']

  
  ########################
  # RUN PROSAIL IN FORWARD MODE
   
  if soil_path is None:
    # Call RTM and generate LUT
    try:
        generate_spectra(
            output_dir=out_dir,
            lut_params=lut_params,
            lut_config=lut_config,
            rtm_config=rtm_config,
            traits=traits
        )

    except Exception as e:
        print(f'Error: {e}')
        pass


  if soil_path is not None:
    # Loop over soil spectra
    soil_df = pd.read_pickle(soil_path)

    try:
        generate_spectra_soil(
            output_dir=out_dir,
            lut_params=lut_params,
            lut_config=lut_config,
            rtm_config=rtm_config,
            traits=traits,
            soil_df=soil_df
        )
    except Exception as e:
        print(f'Error: {e}')
        pass

