import pandas as pd
from ..tools import  toflux

def make_mock(table, sky_center, pixel_scale=0.2, image_offset=(2499.5, 2499.5)):
    """
    Create a galaxy catalog from the output of the KDE resampler

    This is what can be passed to the Frame class

    Parameters
    ----------
    table: pd.Dataframe
        galaxy catalog drawn from the outputs of random samples
    sky_center: tuple
        sky center in ra,dec
    pixel_scale: float
        pixel scale in pixel / arcsec
    image_offset; tuple
        Center of the image, should be at the middle of the canvas (5k by 5k canvases fit a cluster nicely)
    Returns
    -------
        Mock catalog which can be passed to the Frame object

    """
    mock_catalog = pd.DataFrame()
    mock_catalog["RA"] = table["ra"]
    mock_catalog["DEC"] = table["dec"]
    mock_catalog["X"] = (table["ra"] - sky_center[0]) * 60 * 60 / pixel_scale + image_offset[0]
    mock_catalog["Y"] = (table["dec"] - sky_center[1]) * 60 * 60 / pixel_scale + image_offset[1]

    # This is where the translation begins
    mock_catalog["MAG_G"] = table["mag_g"]
    mock_catalog["FLUX_G"] = toflux(table["mag_g"])
    mock_catalog["MAG_R"] = table["mag_i"]
    mock_catalog["FLUX_R"] = toflux(table["mag_r"])
    mock_catalog["MAG_I"] = table["mag_i"]
    mock_catalog["FLUX_I"] = toflux(table["mag_i"])
    mock_catalog["MAG_Z"] = table["mag_z"]
    mock_catalog["FLUX_Z"] = toflux(table["mag_z"])

    mock_catalog["TSIZE"] = table["size_true"]  #
    mock_catalog["FRACDEV"] = 0  # This is just a placeholder, put your own function here if you know how that should look like
    mock_catalog["G1"] = table["ellipticity_1_true"]
    mock_catalog["G2"] = table["ellipticity_2_true"]

    return mock_catalog


def make_transformed_curated_cluster(curated_cluster, sky_center, pixel_scale=0.2, image_offset=(2499.5, 2499.5)):
    """
    Create a galaxy catalog from the output of the KDE resampler,

    **This function exists for historical reasons as a backup in notation**

    This is what can be passed to the Frame class

    Parameters
    ----------
    table: pd.Dataframe
        galaxy catalog drawn from the outputs of random samples
    sky_center: tuple
        sky center in ra,dec
    pixel_scale: float
        pixel scale in pixel / arcsec
    image_offset; tuple
        Center of the image, should be at the middle of the canvas (5k by 5k canvases fit a cluster nicely)
    Returns
    -------
        Mock catalog which can be passed to the Frame object

    """
    mock_cluster = pd.DataFrame()
    mock_cluster["RA"] = curated_cluster["RA"] / 60. + sky_center[0]
    mock_cluster["DEC"] = curated_cluster["DEC"] / 60. + sky_center[1]
    mock_cluster["X"] = (curated_cluster["RA"]) * 60 / pixel_scale + image_offset[0]
    mock_cluster["Y"] = (curated_cluster["DEC"]) * 60 / pixel_scale + image_offset[1]

    mock_cluster["MAG_G"] = curated_cluster["MAG_G"]
    mock_cluster["FLUX_G"] = curated_cluster["FLUX_G"]
    mock_cluster["MAG_R"] = curated_cluster["MAG_R"]
    mock_cluster["FLUX_R"] = curated_cluster["FLUX_R"]
    mock_cluster["MAG_I"] = curated_cluster["MAG_I"]
    mock_cluster["FLUX_I"] = curated_cluster["FLUX_I"]
    mock_cluster["MAG_Z"] = curated_cluster["MAG_Z"]
    mock_cluster["FLUX_Z"] = curated_cluster["FLUX_Z"]

    mock_cluster["TSIZE"] = curated_cluster["TSIZE"]  #
    mock_cluster["FRACDEV"] = curated_cluster[
        "FRACDEV"]  # This is just a placeholder, put your own function here if you know how that should look like
    mock_cluster["G1"] = curated_cluster["G1"]
    mock_cluster["G2"] = curated_cluster["G2"]

    return mock_cluster