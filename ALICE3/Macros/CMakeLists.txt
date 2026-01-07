# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

install(FILES Configuration/a3geo.ini
              Configuration/a3geometry_v2_10kG.ini
              Configuration/a3geometry_v2_20kG_dipole.ini
              Configuration/a3geometry_v2_20kG.ini
              Configuration/a3geometry_v3.ini
              Configuration/a3geometry_v4.ini
        PERMISSIONS GROUP_READ GROUP_EXECUTE OWNER_EXECUTE OWNER_WRITE OWNER_READ WORLD_EXECUTE WORLD_READ
        DESTINATION share/alice3/)
