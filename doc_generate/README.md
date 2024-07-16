# Documentation Generator

It is a very simple bash script that automatically generates a DRIVEWISE-styled documentation for your Python and/or C/C++ projects. We tried to keep the dependencies at a minimum in an effort to make it as cross-platform as possible.

## Usage

It is very easy, simply follow these steps.

- Run the script with `./docgen.sh` (or `bash docgen.sh` if you do not want to `chmod +x` it);
- Follow the instructions.

## Requirements

We tried to do a good job at checking if you have the required dependencies, and telling you which ones are missing. However, we may have missed something, so this is a complete list of requirements (hopefully):

- [Python](https://www.python.org/) (for any project since Sphinx needs it);
- [Conda](https://docs.conda.io/en/latest/) (optional for Python projects);
- [Doxygen](https://www.doxygen.nl/index.html) (for C/C++ projects).

## Notes

Some useful information.

- Once you have generated the documentation, it will generate a `build_html.sh` bash script that you can use to build the documentation in `HTML` format. Just run it with `./build_html.sh` (or `bash build_html.sh` if you do not want to `chmod +x` it). However, if you know what you are doing you can obviously still use the Sphinx commands.
- Each path you provide to the script will be absolute and hard coded in the scripts. This means that if you move the script, you will have to re-generate the documentation (but you can always edit the scripts by hand if you know how to do it). This will likely change in the future to allow for more flexibility.
- Please consider the documentation we generate as a template. This means that it is expected that if you want to add something very specific (e.g., a new static page) you will do it by yourself.
- We strongly recommend that you keep all your custom files separate from the ones we automatically generated. This way you will be able to easily run this script again to re-generate the template without losing any of your files by accident.

## To-do

- [ ] Test on Windows (using WSL).
- [ ] Reduce to a minimum the absolute paths to have a more flexible `docs` folder.
- [ ] Work on the documentation template (e.g., add a logo, change the colour scheme, change the index page, etc.).
- [ ] Use regular expressions whenever possible to improve robustness.
- [ ] Triple check every `rm` to ensure it does not accidentaly delete everything it shouldn't.

## Warning

This script may work, by it is far from being perfect. At the moment it was only tested on macOS and Linux with just a couple of projects. If you find any bug, please report it, and we can try to fix it together.

## License

We licensed this project under the BSD-2-Clause License — see the [LICENSE](LICENSE) file for details.
