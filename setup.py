import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ConSReg",
    version="1.1.5",
    author="Qi Song",
    author_email="alexsong@vt.edu",
    description="condition-specific regulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LiLabAtVT/ConSReg",
    license='MIT',
    python_requires='>=3.6.0',
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy==1.16.2',
        'scipy==1.1.0',
        'pandas==0.21.1',
        'joblib>=0.11',
        'rpy2==2.8.6',
        'networkx>=2',
        'scikit-learn==0.19.1',
        'intervaltree==2.1.0'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

