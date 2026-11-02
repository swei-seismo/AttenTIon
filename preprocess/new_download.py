import obspy
from obspy.clients.fdsn.mass_downloader import RectangularDomain, Restrictions, MassDownloader
import logging


def download_data(orid):
    ## download data based on picks.dat
    picks_info = "extra_download_data_list.dat"
    with open(picks_info, "r") as f:
        lines = f.readlines()

    picks = [line.strip() for line in lines if orid in line]

    origin_time = obspy.UTCDateTime(picks[0].split()[2])
    start_time = origin_time-60
    end_time = origin_time+300
    stations = set()
    for pick in picks:
        parts = pick.split()
        sta = parts[1]
        stations.add(sta)

    domain = RectangularDomain(minlatitude=51, maxlatitude=60, minlongitude=-163, maxlongitude=-147)
    mdl = MassDownloader(providers=["IRIS"])
    logger = logging.getLogger("obspy.clients.fdsn.mass_downloader")
    logger.setLevel(logging.WARNING)  

    for sta in sorted(stations):
        restrictions = Restrictions(
            starttime=start_time,
            endtime=end_time,
            network="*",
            station=sta)
    
        mseed_storage = f"./waveforms/{orid}/" + "{network}.{station}.{channel}.{location}.{starttime}.{endtime}.mseed"
        stationxml_storage = "./stations"
        try:
            mdl.download(domain, restrictions, mseed_storage=mseed_storage, stationxml_storage=stationxml_storage)
        except Exception as e:
                print(f"Download failed for {sta}: {e}")


with open("extra_Eventname", "r") as f:
    oridlst = [line.strip() for line in f]

for orid in oridlst:
    download_data(orid)