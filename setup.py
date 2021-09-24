from setuptools import setup

setup(
    name='SimBa2',
    version='1.1',
    description='Predicts protein stability changes upon mutation',
    author='Kristoffer T. Bæk',
    author_email='krisb@kemi.dtu.dk',
    url='https://github.com/ktbaek/simba2/',
    packages=['simba2'],
    install_requires=[
        'pandas>=1.2',
        'freesasa>=2.1.0',
        'biopython>=1.7',
        'click>=7.1.2',
        'natsort>=7.1.1'
    ],
    python_requires='>=3.6',
    package_data={'simba2': ['HVdiff_table.json', 'naccess.config']},
    include_package_data = True,
    entry_points={
        'console_scripts': [
            'simba2=simba2.__main__:main'
        ]
    }
)
