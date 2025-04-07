import numpy as np
from astroquery.ipac.ned import Ned

def find_redshift(src_name: str):
    result_table = Ned.query_object(src_name)
    # Check if there are results
    if result_table is not None:
        # Extract the redshift from the result table
        z = float(result_table['Redshift'].data[0])
        if np.isnan(z):
            raise NotImplementedError(f'Failed to find redshift for {src_name}')
    else:
        raise ValueError(f"Object named {src_name} not found in NED database")
    return z
