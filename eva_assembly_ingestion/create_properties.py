# Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.config_utils import get_properties_from_xml_file, get_primary_mongo_creds_for_profile, \
    get_accession_pg_creds_for_profile


def write_remapping_process_props_template(template_file_path):
    mongo_host, mongo_user, mongo_pass = get_primary_mongo_creds_for_profile(cfg['maven']['environment'],
                                                                             cfg['maven']['settings_file'])
    pg_url, pg_user, pg_pass = get_accession_pg_creds_for_profile(cfg['maven']['environment'],
                                                                  cfg['maven']['settings_file'])
    with open(template_file_path, 'w') as open_file:
        open_file.write(f'''spring.datasource.driver-class-name=org.postgresql.Driver
spring.datasource.url={pg_url}
spring.datasource.username={pg_user}
spring.datasource.password={pg_pass}
spring.datasource.tomcat.max-active=3

spring.jpa.generate-ddl=true

spring.data.mongodb.host={mongo_host}
spring.data.mongodb.port=27017
spring.data.mongodb.database=eva_accession_sharded
spring.data.mongodb.username={mongo_user}
spring.data.mongodb.password={mongo_pass}

spring.data.mongodb.authentication-database=admin
mongodb.read-preference=secondaryPreferred
spring.main.web-application-type=none
spring.main.allow-bean-definition-overriding=true
spring.jpa.properties.hibernate.jdbc.lob.non_contextual_creation=true
logging.level.uk.ac.ebi.eva.accession.remapping=INFO
parameters.chunkSize=1000

''')
    return template_file_path


def write_clustering_props_template(template_file_path, instance):
    # Additional properties needed for clustering; these need to be appended to the above.
    properties = get_properties_from_xml_file(cfg['maven']['environment'], cfg['maven']['settings_file'])
    counts_url = properties['eva.count-stats.url']
    counts_username = properties['eva.count-stats.username']
    counts_password = properties['eva.count-stats.password']
    with open(template_file_path, 'w') as open_file:
        open_file.write(f'''
accessioning.instanceId=instance-{instance}
accessioning.submitted.categoryId=ss
accessioning.clustered.categoryId=rs

accessioning.monotonic.ss.blockSize=100000
accessioning.monotonic.ss.blockStartValue=5000000000
accessioning.monotonic.ss.nextBlockInterval=1000000000
accessioning.monotonic.rs.blockSize=100000
accessioning.monotonic.rs.blockStartValue=3000000000
accessioning.monotonic.rs.nextBlockInterval=1000000000

eva.count-stats.url={counts_url}
eva.count-stats.username={counts_username}
eva.count-stats.password={counts_password}
''')
    return template_file_path
