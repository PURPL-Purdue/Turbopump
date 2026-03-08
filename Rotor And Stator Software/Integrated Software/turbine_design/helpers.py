import pandas as pd


def export_parameters_to_excel(file_name, configurations, parameters):
    """
    Exports parameters to an Excel file in the required format.

    Parameters:
        file_name (str): The name of the Excel file to create.
        configurations (list of str): List of configuration names.
        parameters (dict): Dictionary where keys are parameter names and values are lists of parameter values.
    """
    # Prepare the data for the Excel sheet
    param_names = list(parameters.keys())  # Parameter names
    param_values = pd.DataFrame(parameters)  # Parameter values as a DataFrame

    # Insert configurations as the first column
    param_values.insert(0, 'Configurations', configurations)

    # Create an empty DataFrame for the top two blank rows
    blank_rows = pd.DataFrame([[""] * (len(param_names) + 1)] * 2)

    # Combine blank rows, parameter names (header), and parameter data
    header_row = [[""] + param_names]  # Parameter names start from B2
    full_data = pd.concat(
        [blank_rows, pd.DataFrame(header_row), param_values],
        ignore_index=True
    )

    # Save to an Excel file
    with pd.ExcelWriter(file_name, engine='openpyxl') as writer:
        full_data.to_excel(writer, index=False, header=False)

def import_parameters_from_excel(file_name):
    """
    Imports parameters from an Excel file in the required format.

    Parameters:
        file_name (str): The name of the Excel file to read.

    Returns:
        configurations (list of str): List of configuration names.
        parameters (dict): Dictionary where keys are parameter names and values are lists of parameter values.
    """
    df = pd.read_excel(file_name, header=None)

    # Extract configurations and parameters
    configurations = df.iloc[2:, 0].tolist()
    parameters = {col: df.iloc[2:, i].tolist() for i, col in enumerate(df.iloc[1, 1:], start=1)}

    return configurations, parameters

def export_parameters_to_csv(file_name, configurations, parameters):
    """
    Exports parameters to a CSV file in the required format.

    Parameters:
        file_name (str): The name of the CSV file to create.
        configurations (list of str): List of configuration names.
        parameters (dict): Dictionary where keys are parameter names and values are lists of parameter values.
    """
    # Prepare the data for the CSV file
    param_names = list(parameters.keys())  # Parameter names
    param_values = pd.DataFrame(parameters)  # Parameter values as a DataFrame

    # Insert configurations as the first column
    param_values.insert(0, 'Configurations', configurations)

    # Create an empty DataFrame for the top two blank rows
    blank_rows = pd.DataFrame([[""] * (len(param_names) + 1)] * 2)

    # Combine blank rows, parameter names (header), and parameter data
    header_row = [[""] + param_names]  # Parameter names start from B2
    full_data = pd.concat(
        [blank_rows, pd.DataFrame(header_row), param_values],
        ignore_index=True
    )

    # Save to a CSV file
    full_data.to_csv(file_name, index=False, header=False)

def import_parameters_from_csv(file_name):
    """
    Imports parameters from a CSV file in the required format.

    Parameters:
        file_name (str): The name of the CSV file to read.

    Returns:
        configurations (list of str): List of configuration names.
        parameters (dict): Dictionary where keys are parameter names and values are lists of parameter values.
    """
    df = pd.read_csv(file_name, header=None)

    # Extract configurations and parameters
    configurations = df.iloc[2:, 0].tolist()
    parameters = {col: df.iloc[2:, i].tolist() for i, col in enumerate(df.iloc[1, 1:], start=1)}

    return configurations, parameters