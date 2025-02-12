import os
import yaml
import importlib.resources as pkg_resources


def load_config():
    config = {}

    # 1. Load user-specific config
    user_config_path = os.path.expanduser("~/.config/mutpredpy/config.yaml")
    if os.path.exists(user_config_path):
        with open(user_config_path, "r") as file:
            config.update(yaml.safe_load(file))

    # 2. Load packaged default config (as fallback)
    with pkg_resources.files("mutpredpy").joinpath("config.yaml").open("r") as file:
        default_config = yaml.safe_load(file)
        config = {**default_config, **config}  # Merge defaults with user/system settings

    return config