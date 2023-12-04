# Contributing Guidelines

Welcome to python-sscha! We appreciate your interest in contributing. Below are the guidelines to contribute effectively to our codebase.

## Forking and Pull Requests

1. **Fork the Repository:** Start by forking the repository to your GitHub account. This creates a copy where you can make your changes.
2. **Create a Branch:** Within your fork, create a new branch for your changes. This keeps the work organized and separate from the main codebase.
3. **Open a Pull Request:** Once your changes are complete, open a pull request against the main repository. This is where your contributions will be reviewed before merging.

## Making Atomic Changes

- **Keep Changes Well-Localized:** Ensure that your changes are atomic. This means that each change should be self-contained and impact only a specific aspect of the code.
- **Avoid Complete Rewrites:** Use existing code as a basis for your changes. Avoid rewriting everything or duplicating existing code unless absolutely necessary.
- **Make use of CellConstructor** Most utility functions are implemented in CellConstructor, please, look into the code, specifically the Methods.py, Phonons.py and Structure.py to find helper and utility functions.
- **Mind the final user** Do not expect the user to be an expert on what you are doing. Ideally, the final API for the user should be as much as black box as possible. Avoid making the user to do unit conversion from the standard SSCHA and CellConstructor units. 
- **Ensure backward compatibility** Its very bad if an user updates the code and they cannot reuse all their scripts. Ideally, your change shoud not break any existing feature or require the user to change their scripts.

## Adding Tests

- **Write Complete Tests:** For every edit or addition, include a complete test. This ensures that your changes work as expected and do not break existing functionality.
- **Tests inside the tests/ directory** Look for examples on other tests. Include files to reproduce your tests. The tests must be quick to run, maximum few seconds.
- **Run the testsuit yourself** You can install pytest and run all the tests yourself to ensure everything is working before opening the pull request.

## Documentation and Examples

- **Document Your Changes:** Clearly document any changes or additions you make to the codebase. This helps others understand what you've done and why. Use docstrings for functions! (we use reStructuredText .rst format)
- **Provide Usage Examples:** Include examples demonstrating how to use the new implementations or changes. Put an example inside the Examples/ directory.
- **Write a tutorial** The best way to have your pull request quickly merged into the main code is to write a usefull tutorial for users on how to use your new functionality.

By following these guidelines, you help maintain the quality and integrity of the SSCHA project. We look forward to your contributions!

