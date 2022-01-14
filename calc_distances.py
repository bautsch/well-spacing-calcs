import pandas as pd
import numpy as np
from scipy.stats import mode

survey_path = r'survey_data.csv'
well_path = r'well_dates.csv'

survey_df = pd.read_csv(survey_path, index_col='well_id')
well_df = pd.read_csv(well_path, parse_dates=['compl_date'], index_col='well_id')
well_count = len(well_df)
parameters = {
                'xpos_dist': np.nan,
                'xpos_date': pd.NaT,
                'xpos_id': np.nan,
                'xpos_api': np.nan,
                'xneg_dist': np.nan,
                'xneg_date': pd.NaT,
                'xneg_id': np.nan,
                'xneg_api': np.nan,
                'ypos_dist': np.nan,
                'ypos_date': pd.NaT,
                'ypos_id': np.nan,
                'ypos_api': np.nan,
                'yneg_dist': np.nan,
                'yneg_date': pd.NaT,
                'yneg_id': np.nan,
                'yneg_api': np.nan,
                'zpos_dist': np.nan,
                'zpos_date': pd.NaT,
                'zpos_id': np.nan,
                'zpos_api': np.nan,
                'zneg_dist': np.nan,
                'zneg_date': pd.NaT,
                'zneg_id': np.nan,
                'zneg_api': np.nan,
                'xy_dist': np.nan,
                'xy_date': pd.NaT,
                'xy_id': np.nan,
                'xy_api': np.nan,
                'xyz_dist': np.nan,
                'xyz_date': pd.NaT,
                'xyz_id': np.nan,
                'xyz_api': np.nan
             }

# Filter to lateral only (80 deg of inclination or greater)
lateral_df = survey_df[survey_df['deg_inclination'] >= 88]

xmin = lateral_df.groupby(by=['well_id'])['x'].min()
xmax = lateral_df.groupby(by=['well_id'])['x'].max()
xmid = ((xmin + xmax) / 2.).rename('xmid', inplace=True)
xlen = (xmax - xmin).rename('xlen', inplace=True)

ymin = lateral_df.groupby(by=['well_id'])['y'].min()
ymax = lateral_df.groupby(by=['well_id'])['y'].max()
ymid = ((ymin + ymax) / 2.).rename('ymid', inplace=True)
ylen = (ymax - ymin).rename('ylen', inplace=True)

zmin = lateral_df.groupby(by=['well_id'])['z'].min()
zmax = lateral_df.groupby(by=['well_id'])['z'].max()
zmid = ((zmin + zmax) / 2.).rename('zmid', inplace=True)
zlen = (zmax - zmin).rename('zlen', inplace=True)

well_df = well_df.join([xmid, xlen, ymid, ylen, zmid, zlen])

well_df['NS_orientation'] = ylen > xlen

# Generate list of sorted, unique first production dates to use as time steps
date_list = well_df['compl_date'].unique()
date_list = np.sort(date_list)
i = 0
wells = set()
for date in date_list:
    # if i > 120:
    #     break
    print('date', date)
    survey_subset = well_df[well_df['compl_date'] <= date]
    survey_subset = survey_subset[~pd.isna(survey_subset.xmid)]
    well_ids = survey_subset.index
    for well in well_ids:
        if well in wells:
            continue
        if pd.isna(well_df.loc[well].xmid):
            for param in parameters:
                if 'date' in param:
                    well_df.loc[well, param] = pd.NaT
                else:
                    well_df.loc[well, param] = np.nan        
            continue
        elif survey_subset[survey_subset.index != well].xmid.isnull().values.all():
            for param in parameters:
                if 'date' in param:
                    well_df.loc[well, param] = pd.NaT
                else:
                    well_df.loc[well, param] = np.nan 
            continue
        else:
            distance_df = survey_subset[survey_subset.index != well].copy()
            distance_df.loc[:, 'xmid'] -= well_df.loc[well, 'xmid']
            distance_df.loc[:, 'ymid'] -= well_df.loc[well, 'ymid']
            distance_df.loc[:, 'zmid'] -= well_df.loc[well, 'zmid']
            distance_df.loc[:, 'xymid'] = np.sqrt(np.power(distance_df['xmid'], 2)
                                            + np.power(distance_df['ymid'], 2))
            distance_df.loc[:, 'xyzmid'] = np.sqrt(np.power(distance_df['xmid'], 2)
                                            + np.power(distance_df['ymid'], 2)
                                            + np.power(distance_df['zmid'], 2))

            if well_df.loc[well, 'NS_orientation']:
                xpos = distance_df[(distance_df.xmid>0) 
                                   & (abs(distance_df.ymid)<distance_df.ylen)
                                   & (abs(distance_df.zmid)<300)].xmid.min()
                xneg = distance_df[(distance_df.xmid<0) 
                                   & (abs(distance_df.ymid)<distance_df.ylen)
                                   & (abs(distance_df.zmid)<300)].xmid.max()
                ypos = distance_df[(distance_df.ymid>0) 
                                   & (abs(distance_df.xmid)<1000)
                                   & (abs(distance_df.zmid)<300)].ymid.min()
                yneg = distance_df[(distance_df.ymid<0) 
                                   & (abs(distance_df.xmid)<1000)
                                   & (abs(distance_df.zmid)<300)].ymid.max()
            else:
                xpos = distance_df[(distance_df.xmid>0)
                                   & (abs(distance_df.ymid)<1000)
                                   & (abs(distance_df.zmid)<300)].xmid.min()
                xneg = distance_df[(distance_df.xmid<0) 
                                   & (abs(distance_df.ymid)<1000)
                                   & (abs(distance_df.zmid)<300)].xmid.max()
                ypos = distance_df[(distance_df.ymid>0) 
                                   & (abs(distance_df.xmid)<distance_df.xlen)
                                   & (abs(distance_df.zmid)<300)].ymid.min()
                yneg = distance_df[(distance_df.ymid<0) 
                                   & (abs(distance_df.xmid)<distance_df.xlen)
                                   & (abs(distance_df.zmid)<300)].ymid.max()
                    
            zpos = distance_df[(distance_df.zmid>0) & (abs(distance_df.xymid)<1000)].zmid.min()
            zneg = distance_df[(distance_df.zmid<0) & (abs(distance_df.xymid)<1000)].zmid.max()
            xy = distance_df.xymid.min()
            xyz = distance_df.xyzmid.min()

            xy_idx = distance_df[distance_df.xymid == xy].index[0]
            xy_date = distance_df.loc[xy_idx].compl_date
            xy_api = distance_df.loc[xy_idx].api_formatted
            xyz_idx = distance_df[distance_df.xyzmid == xyz].index[0]
            xyz_date = distance_df.loc[xyz_idx].compl_date
            xyz_api = distance_df.loc[xyz_idx].api_formatted

            if pd.isna(xpos):
                xpos_date = pd.NaT
                xpos_idx = np.nan
                xpos_api = np.nan
            else:
                xpos_idx = distance_df[distance_df.xmid == xpos].index[0]
                xpos_date = distance_df.loc[xpos_idx].compl_date
                xpos_api = distance_df.loc[xpos_idx].api_formatted

            if pd.isna(xneg):
                xneg_date = pd.NaT
                xneg_idx = np.nan
                xneg_api = np.nan
            else:
                xneg_idx = distance_df[distance_df.xmid == xneg].index[0]
                xneg_date = distance_df.loc[xneg_idx].compl_date
                xneg_api = distance_df.loc[xneg_idx].api_formatted

            if pd.isna(ypos):
                ypos_date = pd.NaT
                ypos_idx = np.nan
                ypos_api = np.nan
            else:
                ypos_idx = distance_df[distance_df.ymid == ypos].index[0]
                ypos_date = distance_df.loc[ypos_idx].compl_date
                ypos_api = distance_df.loc[ypos_idx].api_formatted

            if pd.isna(yneg):
                yneg_date = pd.NaT
                yneg_idx = np.nan
                yneg_api = np.nan
            else:
                yneg_idx = distance_df[distance_df.ymid == yneg].index[0]
                yneg_date = distance_df.loc[yneg_idx].compl_date
                yneg_api = distance_df.loc[yneg_idx].api_formatted

            if pd.isna(zpos):
                zpos_date = pd.NaT
                zpos_idx = np.nan
                zpos_api = np.nan
            else:
                zpos_idx = distance_df[distance_df.zmid == zpos].index[0]
                zpos_date = distance_df.loc[zpos_idx].compl_date
                zpos_api = distance_df.loc[zpos_idx].api_formatted

            if pd.isna(zneg):
                zneg_date = pd.NaT
                zneg_idx = np.nan
                zneg_api = np.nan
            else:
                zneg_idx = distance_df[distance_df.zmid == zneg].index[0]
                zneg_date = distance_df.loc[zneg_idx].compl_date
                zneg_api = distance_df.loc[zneg_idx].api_formatted

            parameters['xpos_dist'] = xpos
            parameters['xpos_date'] = xpos_date
            parameters['xpos_id'] = xpos_idx
            parameters['xpos_api'] = xpos_api
            parameters['xneg_dist'] = xneg
            parameters['xneg_date'] = xneg_date
            parameters['xneg_id'] = xneg_idx
            parameters['xneg_api'] = xneg_api
            parameters['ypos_dist'] = ypos
            parameters['ypos_date'] = ypos_date
            parameters['ypos_id'] = ypos_idx
            parameters['ypos_api'] = ypos_api
            parameters['yneg_dist'] = yneg
            parameters['yneg_date'] = yneg_date
            parameters['yneg_id'] = yneg_idx
            parameters['yneg_api'] = yneg_api
            parameters['zpos_dist'] = zpos
            parameters['zpos_date'] = zpos_date
            parameters['zpos_id'] = zpos_idx
            parameters['zpos_api'] = zpos_api
            parameters['zneg_dist'] = zneg
            parameters['zneg_date'] = zneg_date
            parameters['zneg_id'] = zneg_idx
            parameters['zneg_api'] = zneg_api
            parameters['xy_dist'] = xy
            parameters['xy_date'] = xy_date
            parameters['xy_api'] = xy_api
            parameters['xy_id'] = xy_idx
            parameters['xyz_dist'] = xyz
            parameters['xyz_api'] = xyz_api
            parameters['xyz_date'] = xyz_date
            parameters['xyz_id'] = xyz_idx

            for param in parameters:
                well_df.loc[well, param] = parameters[param]

    wells.update(list(well_ids))
    i += 1

well_df.to_csv('well_distances/2018-11-20/well_distances.csv')