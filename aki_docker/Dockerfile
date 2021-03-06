FROM dbmi/4ce-analysis

## customization for 4CE Phase 2.0 specific for AKI project
RUN R -e "install.packages('data.table')"
RUN R -e "install.packages('tidyverse')"
RUN R -e "install.packages('plyr')"
RUN R -e "install.packages('exactmeta')"
RUN R -e "install.packages('fmsb')"
RUN R -e "install.packages('forestplot')"
RUN R -e "install.packages('httr')"
RUN R -e "install.packages('lme4')"
RUN R -e "install.packages('metafor')"
RUN R -e "install.packages('rtf')"
RUN R -e "install.packages('splines')"
RUN R -e "install.packages('survival')"
RUN R -e "install.packages('RcppRoll')"
RUN R -e "install.packages('zoo')"
RUN R -e "install.packages('survminer')"

## install RStudio Server / notebooks
WORKDIR /opt/rstudioserver

RUN wget http://archive.ubuntu.com/ubuntu/pool/main/o/openssl1.0/libssl1.0.0_1.0.2n-1ubuntu5.4_amd64.deb
RUN dpkg -i ./libssl1.0.0_1.0.2n-1ubuntu5.4_amd64.deb

## RUN apt-get update
RUN apt-get install -y gdebi-core
RUN wget https://download2.rstudio.org/server/trusty/amd64/rstudio-server-1.2.5042-amd64.deb
RUN gdebi -n rstudio-server-1.2.5042-amd64.deb

## Copy startup script.
COPY startup.sh /startup/startup.sh
RUN chmod 700 /startup/startup.sh

WORKDIR /4ceData

## Lock the default user from analysis docker
RUN usermod -L dockeruser
RUN chage -E0 dockeruser
RUN usermod -s /sbin/nologin dockeruser

# Copy RStudio Config
COPY rserver.conf /etc/rstudio/rserver.conf

CMD ["/startup/startup.sh"]