# Experiment 0

PH distribution estimation using EM algorithm


## Build

```bash
docker build -t empht .
```

## Usage

```bash
sh ./exec.sh <phases> <data_file> <iterations>
```

### Parameters

- `phases`: Number of phases for the PH distribution
- `data_file`: Path to the input data file
- `iterations`: Number of EM algorithm iterations

### Example

```bash
sh ./exec.sh 10 ./data/unweighted100 10
```

