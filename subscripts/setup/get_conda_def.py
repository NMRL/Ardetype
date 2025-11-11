import subprocess
import yaml
import json
import argparse

def fetch_package_dependencies(package_name, package_version, channels):
    """Fetch dependencies of a given package using mamba search."""
    try:
        # Construct the mamba search command
        command = ["mamba", "search", "--json", f"{package_name}={package_version}"]
        for channel in channels:
            command.extend(["-c", channel])
        
        # Execute the mamba search command
        result = subprocess.run(command, stdout=subprocess.PIPE, check=True)
        
        # Parse the JSON output
        package_info = json.loads(result.stdout)
        package_info = package_info['result']['pkgs'][0]

        # Extract dependencies
        if package_name in package_info['fn']:
            dependencies = package_info.get("depends", [])
        else:
            dependencies = []

        return dependencies

    except subprocess.CalledProcessError as e:
        print(f"Error while fetching dependencies: {e}")
        return []

def generate_environment_yaml(package_name, version, dependencies, channels, output_file="environment.yml"):
    """Generate a Conda environment YAML file."""
    env_data = {
        "name": f"{package_name}_{version}",
        "channels": channels,
        "dependencies": [f"{package_name}={version}"] + dependencies,
    }

    with open(output_file, "w") as yaml_file:
        yaml.dump(env_data, yaml_file, default_flow_style=False)

    print(f"Environment YAML generated: {output_file}")
    return output_file

def main():
    """Main function to fetch metadata and generate YAML."""
    parser = argparse.ArgumentParser(description="Generate a Conda environment YAML file for a package.")
    parser.add_argument("--package", required=True, help="Name of the package")
    parser.add_argument("--version", required=True, help="Version of the package")
    parser.add_argument("--output", required=True, help="Output YAML file path")
    parser.add_argument("--channels", nargs='+', required=False, help="List of channels to search", default=["bioconda","conda-forge","defaults"])
    

    args = parser.parse_args()

    package_name = args.package
    version = args.version
    channels = args.channels
    output_file = args.output

    print(f"Fetching dependencies for {package_name} from {channels}...")
    dependencies = fetch_package_dependencies(package_name=package_name, package_version=version, channels=channels)

    if not dependencies:
        print(f"No dependencies found for {package_name}.")
        return

    # Generate environment YAML
    generate_environment_yaml(package_name, version, dependencies, channels, output_file)

if __name__ == "__main__":
    main()