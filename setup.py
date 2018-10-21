import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ConSReg",
    version="1.0.0",
    author="Qi Song",
    author_email="alexsong@vt.edu",
    description="condition-specific regulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LiLabAtVT/ConSReg",
    license='MIT',
    python_requires='>=2.7',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

