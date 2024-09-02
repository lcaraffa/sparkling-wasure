from setuptools import setup

setup(name="convertcloud",
      version="0.1.0",
      description="Pointcloud format converter",
      author="Alvaro Capellan, Olivier Roulet-Dubonnet",
      author_email="capellan.alvaro@gmail.com",
      url="https://github.com/alvcap/convertcloud",
      packages=["convertcloud"],
      provides=["convertcloud"],
      license="GNU Lesser General Public License v3",

      entry_points={"console_scripts":
          ["cvc = convertcloud.converter:main"]
                    },
      classifiers=[
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering"
      ]
      )
