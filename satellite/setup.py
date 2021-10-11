import setuptools

setuptools.setup( name='satellite',
	version='1.3',
	description='Spectroscopic Analysis Tool for intEgraL FieLd unIt daTacubEs',
	url='',
	author='Stavros Akras',
	author_email='stavrosakras@gmail.com',
	packages=setuptools.find_packages(),
	install_requires=['matplotlib','numpy','seaborn','astropy','scipy','PyNeb']
	)	
