"""
Configuration Module for MutPredPy

This module handles loading, parsing, and managing configuration settings 
for MutPredPy, including paths, system parameters, and runtime options.
"""

import importlib.resources as pkg_resources
import os
import yaml


def load_config():
    """
    Loads and parses the configuration settings from a YAML file that is
    either a user configuration file or a package specific config YAML.

    Returns:
        dict: Parsed configuration settings.

    Raises:
        FileNotFoundError: If the specified configuration file is not found.
        yaml.YAMLError: If there is an error in parsing the YAML file.
    """
    config = {}

    # 1. Load user-specific config
    user_config_path = os.path.expanduser("~/.config/mutpredpy/config.yaml")
    if os.path.exists(user_config_path):
        with open(user_config_path, "r", encoding="utf-8") as file:
            config.update(yaml.safe_load(file))

    # 2. Load packaged default config (as fallback)
    with pkg_resources.files("mutpredpy").joinpath("config.yaml").open("r") as file:
        default_config = yaml.safe_load(file)
        config = {
            **default_config,
            **config,
        }  # Merge defaults with user/system settings

    return config
