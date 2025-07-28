# Strainify

Strainify is an accurate strain-level abundance analysis tool for short-read metagenomics.


## Installation


### Steps

```bash
# Clone the repository
git clone https://github.com/treangenlab/Strainify.git
cd Strainify

# Install dependencies
conda env create -f environment.yaml

```


## Usage

### CLI

```bash
snakemake --cores 12 --configfile config.yaml
```


## 📊 Examples

Include example usage, sample input/output, or screenshots here.

```bash
# Example terminal output
python main.py --demo
```

```text
Output:
Hello from Project Name!
```

---

## Configuration

You can configure the app using CLI flags, environment variables, or config files.

```bash
python main.py --config config.yaml
```

```yaml
# Example config.yaml
param1: value1
param2: value2
```


## Testing

Run the test with:

```bash
test example/single
```

Or:

```bash
test example/paired
```


## Questions / Contact

For questions or suggestions, open an issue or contact [rl152@rice.edu](mailto:rl152@rice.edu).

