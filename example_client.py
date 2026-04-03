#!/usr/bin/env python3
"""
Example client for the TADF Rate Constants Calculator API.
This script demonstrates how to call the various endpoints of the API.
"""

import requests
import json
import matplotlib.pyplot as plt
import numpy as np
import os

# Base URL of the API
BASE_URL = "http://localhost:8000"

def print_header(text):
    """Helper function to print a header"""
    print("\n" + "="*50)
    print(text)
    print("="*50)

def main():
    # Example 1: Convert lifetime to rate constant
    print_header("Example 1: Convert lifetime to rate constant")
    tau_values = [10e-9, 2.5e-6]  # 10 ns, 2.5 μs
    
    response = requests.post(
        f"{BASE_URL}/tau2k",
        json={"tau": tau_values}
    )
    
    if response.status_code == 200:
        result = response.json()
        print(f"Lifetimes: {tau_values} seconds")
        print(f"Rate constants: {result['k']} s^-1")
    else:
        print(f"Error: {response.status_code}")
        print(response.text)

    # Example 2: Calculate phi_PF and phi_DF
    print_header("Example 2: Calculate phi_PF and phi_DF")
    
    data = {
        "PLQY": 0.82,
        "tauPF": 15e-9,  # 15 ns
        "tauDF": 2.9e-6,  # 2.9 μs
        "B_PF": 0.7,
        "B_DF": 0.3
    }
    
    response = requests.post(
        f"{BASE_URL}/cal_phi_PF_DF",
        json=data
    )
    
    if response.status_code == 200:
        result = response.json()
        print(f"PLQY: {data['PLQY']}")
        print(f"Phi_PF: {result['phi_PF']}")
        print(f"Phi_DF: {result['phi_DF']}")
    else:
        print(f"Error: {response.status_code}")
        print(response.text)

    # Example 3: Calculate determined intrinsic rate constants
    print_header("Example 3: Calculate determined intrinsic rate constants")
    
    kPF = 1/data["tauPF"]  # Convert lifetime to rate constant
    kDF = 1/data["tauDF"]  # Convert lifetime to rate constant
    
    # Get phi_PF and phi_DF from Example 2
    phi_PF = result['phi_PF']
    phi_DF = result['phi_DF']
    
    data = {
        "kPF": kPF,
        "kDF": kDF,
        "phi_PF": phi_PF,
        "phi_DF": phi_DF
    }
    
    response = requests.post(
        f"{BASE_URL}/cal_determined_intrinsic_rate_constants",
        json=data
    )
    
    if response.status_code == 200:
        result = response.json()
        print(f"ksr: {result['ksr']:.4e} s^-1")
        print(f"ks: {result['ks']:.4e} s^-1")
        print(f"kt: {result['kt']:.4e} s^-1")
        print(f"kisckrisc: {result['kisckrisc']:.4e} s^-2")
    else:
        print(f"Error: {response.status_code}")
        print(response.text)

    # Example 4: Calculate intrinsic rate constants
    print_header("Example 4: Calculate intrinsic rate constants")
    
    # Create an array of phi_Tnr_PL values ranging from 0 to 1-PLQY
    phi_Tnr_PL_values = np.linspace(0, 1-data["PLQY"], 5).tolist()  # Using 5 points for simplicity
    
    data_intrinsic = {
        "kPF": kPF,
        "kDF": kDF,
        "phi_PF": phi_PF,
        "phi_DF": phi_DF,
        "phi_Tnr_PL": phi_Tnr_PL_values
    }
    
    response = requests.post(
        f"{BASE_URL}/cal_intrinsic_rate_constants",
        json=data_intrinsic
    )
    
    if response.status_code == 200:
        result = response.json()
        print(f"phi_Tnr_PL values: {phi_Tnr_PL_values}")
        print(f"kisc values: {result['kisc']}")
        print(f"krisc values: {result['krisc']}")
    else:
        print(f"Error: {response.status_code}")
        print(response.text)

    # Example 5: Run script
    print_header("Example 5: Run script")
    
    script_data = {
        "tau_PF": data["tauPF"],
        "tau_DF": data["tauDF"],
        "phi_PF": phi_PF,
        "phi_DF": phi_DF,
        "fname": "test_api_results"
    }
    
    response = requests.post(
        f"{BASE_URL}/script",
        json=script_data
    )
    
    if response.status_code == 200:
        result = response.json()
        print("Script executed successfully!")
        print("Results summary:")
        for key, value in result.items():
            if isinstance(value, (int, float)):
                print(f"{key}: {value:.4e}")
            else:
                print(f"{key}: {value}")
    else:
        print(f"Error: {response.status_code}")
        print(response.text)

if __name__ == "__main__":
    main() 