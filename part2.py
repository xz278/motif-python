if __name__ == '__main__':
	import pandas as pd
	import geohash


	print('part2 starts.')

	l = (pd.Series.from_csv('temp2.csv')).tolist()
	x = [geohash.decode(t) for t in l]
	x = pd.DataFrame(x)
	x.to_csv('temp3.csv')

	print('part2 ends.')