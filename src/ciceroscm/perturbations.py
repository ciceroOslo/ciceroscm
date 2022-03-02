"""
Perturbation related methods
"""
import pandas as pd


def perturb_emissions(perturbation_file, emissions_df):
    """
    Adding emissions perturbations to emissions_df
    """
    pert_df = pd.read_csv(perturbation_file, index_col=None)
    for row in pert_df.itertuples(index=True, name="Pandas"):
        print(emissions_df[row.component][row.year])
        emissions_df[row.component][row.year] = (
            emissions_df[row.component][row.year] + row.emission
        )
        print(emissions_df[row.component][row.year])
        print(row.emission)
    return emissions_df
