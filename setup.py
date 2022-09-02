from setuptools import setup

setup (
    entry_points={
        'console_scripts': [
            'streammd = streammd.markdups:main',
            'memcalc = streammd.memcalc:main'
        ]
    }
)
