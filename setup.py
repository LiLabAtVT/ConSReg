import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ConSReg",
    version="1.0.7",
    author="Qi Song",
    author_email="alexsong@vt.edu",
    description="condition-specific regulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LiLabAtVT/ConSReg",
    license='MIT',
    python_requires='>=2.7',
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy>=1.9.0',
        'scipy',
        'pandas==0.21.1',
        'joblib',
        'rpy2==2.8.6',
        'networkx>=2',
        'sklearn',
        'intervaltree'
    ],
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

