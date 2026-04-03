# TADF Rate Constants Calculator API

This is a FastAPI-based backend for the TADF (Thermally Activated Delayed Fluorescence) rate constants calculator. It provides API endpoints to perform various calculations related to exciton dynamics in TADF materials.

## Installation

1. Clone this repository
2. Install the dependencies:
   ```
   pip install -r requirements.txt
   ```

## Running the API

Start the API server with:

```bash
uvicorn app:app --reload
```

The API will be available at http://localhost:8000

## Interactive API Documentation

FastAPI automatically generates interactive API documentation:
- Swagger UI: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc

## API Endpoints

The API provides the following endpoints:

### Basic Conversions

- `POST /tau2k`: Convert lifetime to rate constant
- `POST /k2tau`: Convert rate constant to lifetime

### Quantum Yield Calculations

- `POST /cal_phi_PF_DF`: Calculate the absolute ratio of prompt and delayed fluorescence

### Rate Constants Calculations

- `POST /cal_determined_intrinsic_rate_constants`: Calculate determinable rate constants from experimental data
- `POST /cal_intrinsic_rate_constants`: Calculate all intrinsic rate constants

### Full Script

- `POST /script`: Run the full rate constants calculator script with provided parameters

## Example Usage

Here's how to use the API with Python requests:

```python
import requests
import json

# Example: Convert lifetime to rate constant
response = requests.post(
    "http://localhost:8000/tau2k",
    json={"tau": [10.0, 20.0]}
)
print(response.json())  # Output: {"k": [0.1, 0.05]}

# Example: Calculate phi_PF and phi_DF
response = requests.post(
    "http://localhost:8000/cal_phi_PF_DF",
    json={
        "PLQY": 0.8,
        "tauPF": 10.0,
        "tauDF": 1000.0,
        "B_PF": 0.6,
        "B_DF": 0.4
    }
)
print(response.json())  # Output: {"phi_PF": value, "phi_DF": value}
```

## Using with Frontend Applications

This API can be consumed by any frontend application including:
- Web applications using JavaScript/React/Vue/Angular
- Mobile applications
- Desktop applications
- Other Python scripts or applications

CORS is enabled, allowing the API to be called from browser-based applications on different domains.

