import re
import requests
import user_agent


def query_tns_main_name(sourcename):
    url = f"https://www.wis-tns.org/search?&discovered_period_value=10&discovered_period_units=years&unclassified_at=0&classified_sne=0&include_frb=0&name={sourcename}&name_like=1&isTNS_AT=all&public=all&ra=&decl=&radius=&coords_unit=arcsec&reporting_groupid%5B%5D=null&groupid%5B%5D=null&classifier_groupid%5B%5D=null&objtype%5B%5D=null&at_type%5B%5D=null&date_start%5Bdate%5D=&date_end%5Bdate%5D=&discovery_mag_min=&discovery_mag_max=&internal_name=&discoverer=&classifier=&spectra_count=&redshift_min=&redshift_max=&hostname=&ext_catid=&ra_range_min=&ra_range_max=&decl_range_min=&decl_range_max=&discovery_instrument%5B%5D=null&classification_instrument%5B%5D=null&associated_groups%5B%5D=null&official_discovery=0&official_classification=0&at_rep_remarks=&class_rep_remarks=&frb_repeat=all&frb_repeater_of_objid=&frb_measured_redshift=0&frb_dm_range_min=&frb_dm_range_max=&frb_rm_range_min=&frb_rm_range_max=&frb_snr_range_min=&frb_snr_range_max=&frb_flux_range_min=&frb_flux_range_max=&num_page=50&display%5Bredshift%5D=1&display%5Bhostname%5D=1&display%5Bhost_redshift%5D=1&display%5Bsource_group_name%5D=1&display%5Bclassifying_source_group_name%5D=1&display%5Bdiscovering_instrument_name%5D=0&display%5Bclassifing_instrument_name%5D=0&display%5Bprograms_name%5D=0&display%5Binternal_name%5D=1&display%5BisTNS_AT%5D=0&display%5Bpublic%5D=1&display%5Bend_pop_period%5D=0&display%5Bspectra_count%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&display%5Bdiscoverer%5D=1&display%5Bremarks%5D=0&display%5Bsources%5D=0&display%5Bbibcode%5D=0&display%5Bext_catalogs%5D=0"
    agent = user_agent.generate_user_agent()
    headers = {'User-Agent': agent}

    resp = requests.get(url, headers=headers)
    return resp.text


def query_tns_survey_name(sourcename):
    url = f"https://www.wis-tns.org/search?&discovered_period_value=10&discovered_period_units=years&unclassified_at=0&classified_sne=0&include_frb=0&name=&name_like=1&isTNS_AT=all&public=all&ra=&decl=&radius=&coords_unit=arcsec&reporting_groupid[]=null&groupid[]=null&classifier_groupid[]=null&objtype[]=null&at_type[]=null&date_start[date]=&date_end[date]=&discovery_mag_min=&discovery_mag_max=&internal_name={sourcename}&discoverer=&classifier=&spectra_count=&redshift_min=&redshift_max=&hostname=&ext_catid=&ra_range_min=&ra_range_max=&decl_range_min=&decl_range_max=&discovery_instrument[]=null&classification_instrument[]=null&associated_groups[]=null&official_discovery=0&official_classification=0&at_rep_remarks=&class_rep_remarks=&frb_repeat=all&frb_repeater_of_objid=&frb_measured_redshift=0&frb_dm_range_min=&frb_dm_range_max=&frb_rm_range_min=&frb_rm_range_max=&frb_snr_range_min=&frb_snr_range_max=&frb_flux_range_min=&frb_flux_range_max=&num_page=50&display[redshift]=1&display[hostname]=1&display[host_redshift]=1&display[source_group_name]=1&display[classifying_source_group_name]=1&display[discovering_instrument_name]=0&display[classifing_instrument_name]=0&display[programs_name]=0&display[internal_name]=1&display[isTNS_AT]=0&display[public]=1&display[end_pop_period]=0&display[spectra_count]=1&display[discoverymag]=1&display[discmagfilter]=1&display[discoverydate]=1&display[discoverer]=1&display[remarks]=0&display[sources]=0&display[bibcode]=0&display[ext_catalogs]=0&format=html&edit[type]=&edit[objname]=&edit[id]=&sort=asc&order=name"
    agent = user_agent.generate_user_agent()
    headers = {'User-Agent': agent}

    resp = requests.get(url, headers=headers)
    return resp.text


def parse_data(my_text):
    # 
    m = re.findall('(?<=<td class="cell-name"><a href="/object/)(.*)(?=</a></td><td class="cell-reps">)' , my_text)
    if len(m) == 0:
        main_id = None
    else:
        main_id = m[0].split("\">")[1]
        # main_id = m
    

    m = re.findall('(?<=<td class="cell-internal_name">)(.*)(?=</td><td class="cell-groups">)', my_text)
    if len(m) == 0:
        id_ = None
    else:
        id_ = m

    m = re.findall('(?<=<td class="cell-type">)(.*)(?=</td><td class="cell-redshift">)', my_text)
    if len(m) == 0:
        otype = None
    else:
        otype = m

    m = re.search('(?<=<td class="cell-ra">)(.*)(?=</td><td class="cell-decl">)', my_text)
    if m:
        ra = m.group(0)
    else:
        ra = None
     
    m = re.search('(?<=<td class="cell-decl">)(.*)(?=</td><td class="cell-discovery_date">)', my_text)
    if m:
        dec = m.group(0)
    else:
        dec = None
    
    m = re.search('(?<=<td class="cell-discovery_date">)(.*)(?=</td><td class="cell-flux">)', my_text)
    if m:
        discovery_time = m.group(0)
    else:
        discovery_time = None
        
       
    return main_id, id_, otype, ra, dec, discovery_time


def get_id_otype_tns(sourcename):
    my_text = query_tns_survey_name(sourcename)
    return parse_data(my_text)
