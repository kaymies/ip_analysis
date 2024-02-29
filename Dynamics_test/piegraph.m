y2010 = [347/369*100; 14/369*100; 6/369*100; 2/369*100]';
y2011 = [160/369*100; 94/369*100; 106/369*100; 9/369*100]';
labels = {'Agree','Neutral','Disagree','Missing'};

t = tiledlayout(1,2,'TileSpacing','compact');

% Create pie charts
ax1 = nexttile;
pie(ax1,y2010)
% title('Quantifying impairments and outcomes with measurement is important for patient care.')

ax2 = nexttile;
pie(ax2,y2011)
% title('Existing standardized balance measures meet my practice needs.')

% Create legend
lgd = legend(labels);
lgd.Layout.Tile = 'east';