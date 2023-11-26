import yaml
import os

_scriptDir = os.path.dirname(__file__)
configPath = os.path.join(_scriptDir, '../../config.yml')
with open(configPath, 'r') as __file:
    config = yaml.safe_load(__file)