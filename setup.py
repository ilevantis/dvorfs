from setuptools import setup, find_packages

setup(
    name = 'dvorfs',
    version = '1.0.1',
    packages = find_packages(),
    python_requires='>=3.6',
    install_requires=[
        "numpy",
        "pandas",
    ],
    entry_points = {
        'console_scripts': [
            'dvorfs=dvorfs.dvorfs:main',
            'process-genewise=dvorfs.process_genewise:main',
        ],
    }
)
