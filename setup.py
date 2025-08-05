from setuptools import setup, find_packages

setup(
    name="LUTR",
    version="0.2",
    license="GPL3",
    description="UTR extensions for annotations from protein orthology based gene prediction tools using exons from reference based transcriptome assembly",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=["pandas", "numpy"],
    python_requires=">=3.10",
    entry_points={
        "console_scripts": [
            "LUTR=lutr.main:main",
            "assembly=scripts.assembly:main",
            "noUTR=scripts.noUTR:main"
        ],
    },
)
