function [x,y] = convert_to_xy(lat, lon)

dczone = utmzone(mean(lat,'omitnan'),mean(lon,'omitnan'));
utmstruct = defaultm('utm');
utmstruct.zone = dczone;
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

[x, y] = mfwdtran(utmstruct, lat, lon);
