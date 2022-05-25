from setuptools import setup

setup(name='rrap',
      version='0.2.5',
      description='A metagenomic read recruitment data pipeline',
      url='https://github.com/Kojiconner/rrap_metag_pkg/',
      author='Conner Kojima',
      author_email='cykojima@usc.edu',
      license='MIT',
      packages=['rrap'],
      scripts=['bin/rrap', 'bin/rpkm_heater'],
      install_requires=[
        "joblib==0.16.0",
        "scikit-learn==0.23.1",
        "scipy>=1.5.2",
        "threadpoolctl==2.1.0"
      ],
      python_requires=">=3.7.3",
      package_dir={"": "src"},
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      zip_safe=False)
