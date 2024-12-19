import subprocess

def execute_spacer_placer_command(input_file_path, output_directory_path):
    command = ["spacerplacer", input_file_path, output_directory_path]

    # Run the command
    try:
        subprocess.run(command, check=True)
        print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")
    except FileNotFoundError:
        print("The 'spacerplacer' command was not found. Ensure it is installed and in your PATH.")

