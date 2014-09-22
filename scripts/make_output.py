#!/usr/bin/env python

import sys, os, re
from pylab import *
import yaml
import jinja2
import glob
import scipy.stats
from string import *
import gzip
import locale
import logging
import tarfile
import urllib

def findPeakPlots(peaks, idfile, plotdir):

    # idfile:
    # chr19_10764250_10765250 reg1000002
    # chr17_80245500_80246750 reg1000003

    d = {}

    for l in open(idfile):
        t = l.strip().split()
        d[t[1]] = t[0]

    # peaks
    # chr19   10764651        10764790        reg1000002.p1   16.244  +
    # chr19   18133094        18133186        reg1000004.p2   15.881  +

    z = loadtxt(peaks, usecols=[4])
    ids = loadtxt(peaks, usecols=[3], dtype=str)

    zids = zip(z, ids)

    sorted_zids = sorted(zids, key = lambda k: k[0], reverse=True)

    L = len(zids)

    regions = []

    for i in [0.05, 0.1, 0.20, 0.4, 0.6, 0.9]: #[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        regions.append(d[sorted_zids[int(i*L)][1].split('.')[0]])

    plots = []
    for r in regions:
        plots.append(os.path.join(plotdir.rstrip('\.tar\.gz'), r + '.png'))

    return plots


def parseGenes(gene_ids, gene_info):

    ## human and mouse:
    # gene_ids are space separated
    # gene_info is ##-separated

    # gene_ids:  NM_024333 TSC_hg19_v1_chr19_+_4304594 AK315661 AK295523 AK303606 AY032617 BC003124 AK021750 AF316829 DQ891833 DQ895020
    # gene_info: NM_024333|FSD1|79187|fibronectin type III and SPRY domain containing 1##AK315661|FSD1|79187|fibronectin type III and SPRY domain containing 1##AK295523|FSD1|79187|fibronectin type III and SPRY domain containing 1##AK303606|FSD1|79187|fibronectin type III and SPRY domain containing 1##AY032617|FSD1|79187|fibronectin type III and SPRY domain containing 1##BC003124|FSD1|79187|fibronectin type III and SPRY domain containing 1##AK021750|FSD1|79187|fibronectin type III and SPRY domain containing 1

    gene_id_list = gene_ids.split()
    gene_info_list = gene_info.split('##')

    info_dict = {}
    for i in gene_info_list:
        i_t = i.split('|')
        try:
            info_dict[i_t[0]] = '%s|%s' %(i_t[1], i_t[3])
        except IndexError: #this is the case if i is None
            info_dict[i] = i

    gene_list = []

    for i in gene_id_list:
        try:
            gene_list.append(info_dict[i])
        except KeyError:
            continue

    gene_list = unique(gene_list)

    return ' '.join(gene_list)
        

def getSitesDict(names, l):

    d = {}
    for n in names:
        d[n] = []
    for site in l:
        t = site.split(':')
        d[t[0]].append(t[1]+':'+t[2])

    for k in d:
        d[k] = ', '.join(d[k])

    return d


def create_FMIid_dict():

    # datafiles.csv:
    # /import/.../wgEncodeSydhTfbsHelas3Usf2IggmusRawDataRep1.fastq.gz      USF2    wgEncodeSydhTfbsHelas3Usf2IggmusRawDataRep1     No sample description given     /import/.../ARNT_ARNT2_BHLHB2_MAX_MYC_USF1.p2        fastq
    # /import/.../wgEncodeSydhTfbsHelas3Usf2IggmusRawDataRep2.fastq.gz      USF2    wgEncodeSydhTfbsHelas3Usf2IggmusRawDataRep2     No sample description given     /import/.../ARNT_ARNT2_BHLHB2_MAX_MYC_USF1.p2        fastq
    # /import/.../wgEncodeSydhTfbsHelas3InputIggmusSRR502473rep1.fastq.gz   BG      wgEncodeSydhTfbsHelas3InputIggmusSRR502473rep1  No sample description given     None    fastq
    # /import/.../raw.data.fastq/wgEncodeSydhTfbsHelas3InputIggmusSRR502474rep2.fastq.gz   BG      wgEncodeSydhTfbsHelas3InputIggmusSRR502474rep2  No sample description given     None    fastq

    # exisitingSamples.csv
    # /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/samples/wgEncodeSydhTfbsHelas3Usf2IggmusRawDataRep1/USF2_1-fraglen/out_dir/outfile.gz     USF2    /import/bc2/home/nimwegen/GROUP/WMs/Mammals/CurrentWMs/ARNT_ARNT2_BHLHB2_MAX_MYC_USF1.p2
    # /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/samples/wgEncodeSydhTfbsHelas3Usf2IggmusRawDataRep2/USF2_2-fraglen/out_dir/outfile.gz     USF2    /import/bc2/home/nimwegen/GROUP/WMs/Mammals/CurrentWMs/ARNT_ARNT2_BHLHB2_MAX_MYC_USF1.p2
    # /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/samples/wgEncodeSydhTfbsHelas3InputIggmusSRR502473rep1/BG_1-fraglen/out_dir/outfile.gz    BG      None
    # /import/bc2/home/nimwegen/GROUP/DeepSeqPipeline.Piotr/samples/wgEncodeSydhTfbsHelas3InputIggmusSRR502474rep2/BG_2-fraglen/out_dir/outfile.gz    BG      None

    FMIid_dict = {}

    for i, line in enumerate(open('PipelineFiles/datafiles.csv')):
        if i == 0:
            continue
        t = line.strip().split()
        
        fac = t[1]
        try:
            FMIid_dict[fac].append((t[0], t[2], t[5]))
        except KeyError:
            FMIid_dict[fac] = [(t[0], t[2], t[5])]

    for i, line in enumerate(open('PipelineFiles/existingSamples.csv')):
        if i == 0:
            continue
        t = line.strip().split()

        fac = t[1]
        try:
            FMIid_dict[fac].append((t[0], t[0].split('/')[-4], 'shiftedbed'))
        except KeyError:
            FMIid_dict[fac] = [(t[0], t[0].split('/')[-4], 'shiftedbed')]

    return FMIid_dict


def find_Reps_in_FMI(TFs_samples_dict, download_out, FMIpath, outdir):
    """
    This function is actually just used for testing purposes... Finding stuff in FMI repository when using shifted bed files...
    """

    FMIid_dict = create_FMIid_dict()

    TFs_samples_html_dict = {}

    for TF_i in TFs_samples_dict:
        TFdown_out = os.path.join(download_out, TF_i)
        os.mkdir(TFdown_out)
        print TF_i

        for i, samp in enumerate(FMIid_dict[TF_i]): # samp = [fmi_id, file_type]
            print i, samp

            tout1 = os.path.join(TFdown_out, 'Rep_%i' %(i+1))
            os.system('mkdir \'%s\'' %tout1)

            tout1_fromHTML = os.path.join('../../..', tout1)
            rep_html_dict = {}

            rep_html_dict['name'] = '%s Replicate %i' %(TF_i, i+1)
            rep_html_dict['file_type'] = 'fastq' #samp[2]

            # 6: Do it dependent on input file type
            if samp[2] == 'fastq' or samp[2] == 'fasta':
                wig = glob.glob(os.path.join(FMIpath, 'samples', samp[1], '*-wig/out_dir/*'))[0] #use glob to evaluate wildcards
                bed = glob.glob(os.path.join(FMIpath, 'samples', samp[1], '*-bedweight/out_dir/*'))[0]
                os.system('cp %s %s/reads.wig.gz' %(wig, tout1))
                os.system('cp %s %s/reads.bed.gz' %(bed, tout1))
                rep_html_dict['wig'] = os.path.join(tout1_fromHTML, 'reads.wig.gz')
                rep_html_dict['bed'] = os.path.join(tout1_fromHTML, 'reads.bed.gz')


            if not samp[2] == 'shiftedbed' or True:
                try:
                    report = os.path.join(FMIpath, 'samples', samp[1], '*-pdfreport/pdfreport.pdf')
                    os.system('cp %s \'%s\'' %(report, tout1))
                    rep_html_dict['report_pdf'] = os.path.join(tout1_fromHTML, 'pdfreport.pdf')
                except Exception:
                    pass
                try:
                    fragsize_plot = glob.glob(os.path.join(FMIpath, 'samples', samp[1], '*fraglen/FragmentLength_plot.pdf'))[0]
                    os.system('convert \'%s\' \'%s/FragmentLength_plot.png\'' %(fragsize_plot, tout1))
                    rep_html_dict['fragsize_plot'] = os.path.join(tout1_fromHTML, 'FragmentLength_plot.png')
                except Exception:
                    pass

            # 7:
            if samp[2] == 'fastq':
                qf_log = os.path.join(outdir, TF_i + '_%i-qualityfilter/job.stdout' %(i+1))
                for l in open(qf_log):
                    if l.startswith('Original number of reads:'):
                        t = l.strip().lstrip('Original number of reads:').split()
                        rep_html_dict['num_reads_input'] = int(t[0])
                    elif l.startswith('Reads that passed the quality filter:'):
                        t = l.strip().lstrip('Reads that passed the quality filter:').split()
                        rep_html_dict['num_reads_1QC'] = int(t[0])

                trans_log = os.path.join(outdir, TF_i +'_%i-trans/transform_log' %(i+1))
                for l in open(trans_log):
                    if l.startswith('#                              too short :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_tooshort'] = int(t[4])
                    elif l.startswith('#                            low entropy :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_lowentropy'] = int(t[4])
                    elif l.startswith('#                                 passed :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC'] = int(t[3])

                mappable_log = os.path.join(outdir, TF_i + '_%i-mappable/numberMappableReads_log' %(i+1))
                rep_html_dict['num_reads_mapped'] = int(open(mappable_log).read().strip().split()[1])


            elif samp[2] == 'fasta':
                trans_log = os.path.join(outdir, TF_i +'_%i-trans/transform_log' %(i+1))
                for l in open(trans_log):
                    if l.startswith('#                        input sequences :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_input'] = int(t[4])
                    elif l.startswith('#                              too short :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_tooshort'] = int(t[4])
                    elif l.startswith('#                            low entropy :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_lowentropy'] = int(t[4])
                    elif l.startswith('#                                 passed :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC'] = int(t[3])

                mappable_log = os.path.join(outdir, TF_i + '_%i-mappable/numberMappableReads_log' %(i+1))
                rep_html_dict['num_reads_mapped'] = int(open(mappable_log).read().strip().split()[1])


            elif samp[2] == 'bed' or samp[2] == 'shiftedbed':
                # try:
                #     gop = gzip.open(samp[0])
                #     gop.readline()
                #     gop.close()
                #     gzipped = True
                # except IOError:
                #     gop.close()
                #     gzipped = False

                # read_num = 0
                # if gzipped:
                #     for line in gzip.open(samp[0]):
                #         t = line.strip().split()
                #         read_num += float(t[4])
                # else:
                #     for line in open(samp[0]):
                #         t = line.strip().split()
                #         read_num += float(t[4])

                read_num=10000000000
                rep_html_dict['num_reads_input'] = int(read_num)

            try:
                TFs_samples_html_dict[TF_i].append(rep_html_dict)
            except KeyError:
                TFs_samples_html_dict[TF_i] = [rep_html_dict]


    return TFs_samples_html_dict


def find_Reps_in_outdir(TFs_samples_dict, download_out, outdir):

    locale.setlocale(locale.LC_ALL, 'en_US') # for making long integers readable with the help of commas...

    TFs_samples_html_dict = {}

    for TF_i in TFs_samples_dict:
        print TF_i
        TFdown_out = os.path.join(download_out, TF_i)
        os.mkdir(TFdown_out)

        for i, samp in enumerate(TFs_samples_dict[TF_i]): #samp = [sample_path, file_type]
            print i, samp
            tout1 = os.path.join(TFdown_out, 'Rep_%i' %(i+1))
            os.mkdir(tout1)

            tout1_fromHTML = os.path.join('../', tout1)
            rep_html_dict = {}

            rep_html_dict['name'] = '%s Replicate %i' %(TF_i, i+1)
            rep_html_dict['file_type'] = samp[1]

            # 6: Do it dependent on input file type
            if samp[1] == 'fastq' or samp[1] == 'fasta':
                root_name = os.path.split(samp[0])[1].split('.')[0]
                wig = glob.glob(os.path.join(outdir, TF_i + '_%i-wig/out_dir/*%s*' %(i+1, root_name)))[0] #use glob to evaluate wildcards
                bed = glob.glob(os.path.join(outdir, TF_i + '_%i-bedweight/out_dir/*%s*' %(i+1, root_name)))[0]
                os.system('cp %s %s/reads.wig.gz' %(wig, tout1))
                os.system('cp %s %s/reads.bed.gz' %(bed, tout1))
                rep_html_dict['wig'] = os.path.join(tout1_fromHTML, 'reads.wig.gz')
                rep_html_dict['bed'] = os.path.join(tout1_fromHTML, 'reads.bed.gz')


            if not samp[1] == 'shiftedbed':
                report = os.path.join(outdir, TF_i + '_%i-pdfreport/pdfreport.pdf' %(i+1))
                fragsize_plot = glob.glob(os.path.join(outdir, TF_i + '_%i-fraglen/FragmentLength_plot.pdf' %(i+1)))[0]

                os.system('cp %s %s' %(report, tout1))
                os.system('convert %s %s/FragmentLength_plot.png' %(fragsize_plot, tout1))

                rep_html_dict['report_pdf'] = os.path.join(tout1_fromHTML, 'pdfreport.pdf')
                rep_html_dict['fragsize_plot'] = os.path.join(tout1_fromHTML, 'FragmentLength_plot.png')


            # 7:
            if samp[1] == 'fastq':
                qf_log = os.path.join(outdir, TF_i + '_%i-qualityfilter/job.stdout' %(i+1))
                for l in open(qf_log):
                    if l.startswith('Original number of reads:'):
                        t = l.strip().lstrip('Original number of reads:').split()
                        rep_html_dict['num_reads_input'] = locale.format("%d", int(t[0]), grouping=True)
                    elif l.startswith('Reads that passed the quality filter:'):
                        t = l.strip().lstrip('Reads that passed the quality filter:').split()
                        rep_html_dict['num_reads_1QC'] = locale.format("%d", int(t[0]), grouping=True)

                trans_log = os.path.join(outdir, TF_i +'_%i-trans/transform_log' %(i+1))
                for l in open(trans_log):
                    if l.startswith('#                              too short :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_tooshort'] = locale.format("%d", int(t[4]), grouping=True)
                    elif l.startswith('#                            low entropy :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_lowentropy'] = locale.format("%d", int(t[4]), grouping=True)
                    elif l.startswith('#                                 passed :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC'] = locale.format("%d", int(t[3]), grouping=True)

                mappable_log = os.path.join(outdir, TF_i + '_%i-mappable/numberMappableReads_log' %(i+1))
                rep_html_dict['num_reads_mapped'] = locale.format("%d", int(open(mappable_log).read().strip().split()[1]),  grouping=True)


            elif samp[1] == 'fasta':
                trans_log = os.path.join(outdir, TF_i +'_%i-trans/transform_log' %(i+1))
                for l in open(trans_log):
                    if l.startswith('#                        input sequences :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_input'] = locale.format("%d", int(t[4]), grouping=True)
                    elif l.startswith('#                              too short :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_tooshort'] = locale.format("%d", int(t[4]), grouping=True)
                    elif l.startswith('#                            low entropy :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC_lowentropy'] = locale.format("%d", int(t[4]), grouping=True)
                    elif l.startswith('#                                 passed :'):
                        t = l.strip().split()
                        rep_html_dict['num_reads_2QC'] = locale.format("%d", int(t[3]), grouping=True)

                mappable_log = os.path.join(outdir, TF_i + '_%i-mappable/numberMappableReads_log' %(i+1))
                rep_html_dict['num_reads_mapped'] = locale.format("%d", int(open(mappable_log).read().strip().split()[1]), grouping=True)


            elif samp[1] == 'bed' or samp[1] == 'shiftedbed':
                try:
                    gop = gzip.open(samp[0])
                    gop.readline()
                    gop.close()
                    gzipped = True
                except IOError:
                    gop.close()
                    gzipped = False

                read_num = 0
                if gzipped:
                    for line in gzip.open(samp[0]):
                        t = line.strip().split()
                        read_num += float(t[4])
                else:
                    for line in open(samp[0]):
                        t = line.strip().split()
                        read_num += float(t[4])

                rep_html_dict['num_reads_input'] = locale.format("%d", int(read_num), grouping=True)

            try:
                TFs_samples_html_dict[TF_i].append(rep_html_dict)
            except KeyError:
                TFs_samples_html_dict[TF_i] = [rep_html_dict]


    return TFs_samples_html_dict



def main():

    # remove large OUTPUT_FMI directory
    # print 'remove large OUTPUT_FMI directory'
    # os.system('rm -r OUTPUT_FMI')

    # Directory REPORT will contain two directories: downloads and html
    # downloads will contain large files (wig, bed) for downloading
    # html contains everything else plus html pages...

    pipeline_dir = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]

    collout = 'report' #collected output
    if not os.path.exists(collout):
        os.mkdir(collout)

    # load directory with all output html templates into jinja2:
    temp_env = jinja2.Environment(loader=jinja2.FileSystemLoader('%s' %(os.path.join(pipeline_dir, 'templates'))),
                                  trim_blocks=True)

    # Information to collect all output is found in yaml file.
    # To collect preprocessing output I will read datafiles.csv and existingSamples.csv later.
    # read yaml file:
    configfile = sys.argv[1]
    cf = open(configfile)
    params = yaml.load(cf)
    cf.close()

    outdir = params['OUT_DIR']
    FMIpath = params['FMI_PATH']
    organism = params['GENOME']

    intype_dict = dict([('FASTQ_FILES', 'fastq'), ('FASTA_FILES', 'fasta'), ('BED_FILES', 'bed'), ('SHIFTEDBED_FILES', 'shiftedbed')])

    TFs_samples_dict = {} #keys are factor names, values are file paths plus file type
    for intype in ['IP_FASTQ_FILES', 'IP_FASTA_FILES', 'IP_BED_FILES', 'IP_SHIFTEDBED_FILES']: # important: everything is done in the order: FASTQ, FASTA, BED, SHIFTEDBED
        t_dict = params[intype]
        if t_dict == None:
            continue
        for tf in t_dict:
            for sample in t_dict[tf].split():
                try:
                    TFs_samples_dict[tf].append((sample, intype_dict[intype.lstrip('IP_')]))
                except KeyError:
                    TFs_samples_dict[tf] = [(sample, intype_dict[intype.lstrip('IP_')])]

    for intype in ['BG_FASTQ_FILES', 'BG_FASTA_FILES', 'BG_BED_FILES', 'BG_SHIFTEDBED_FILES']:
        t_list = params[intype]
        if t_list == None:
            continue
        for sample in t_list.split():
            try:
                TFs_samples_dict['BG'].append((sample, intype_dict[intype.lstrip('BG_')]))
            except KeyError:
                TFs_samples_dict['BG'] = [(sample, intype_dict[intype.lstrip('BG_')])]

    print TFs_samples_dict

    # 1. Find Z_hist, revcum, peak_plots and numbers of significant regions and number of fitted peaks
    # 2. Find for each WM: Logo, LL-score, AUC, AUC plot, FOV, FOV plot,  enrichment at binding sites, binding sites plot, %True, report pdf
    # 3. Find for denovo and known: loglik contribution plot, motif correlation plot
    # 4. Find report for FgBg
    # 5. Find final file containing everything
    # 6. Find for each sample: wig, bed, report
    # 7. Find for each sample: qualityfilter log (number of input reads, thrown out reads 1st QC round), trans log (thrown out reads second QC round), number mappable reads (summed bedweight weights, actual number of mapped reads)


    # REPORT directory is organised in html and downloads directory. In html there is IP and in downloads there is IP and BG for wig and bed files. 
    # make a directory where all output html files get stored and also pictures and final peaks
    html_out = collout #os.path.join(collout, 'html')
    if not os.path.exists(html_out):
        os.mkdir(html_out)

    download_out = 'downloads' #os.path.join(collout, 'downloads')
    os.mkdir(download_out)

    TFs_samples_html_dict = find_Reps_in_outdir(TFs_samples_dict, download_out, outdir)
    #TFs_samples_html_dict = find_Reps_in_FMI(TFs_samples_dict, download_out, FMIpath, outdir) #I think this function shouldn't be of any use any longer.


    for TFname in TFs_samples_dict:
        if TFname == 'BG':
            continue

        # for html
        TFhtml = collout #os.path.join(html_out, TFname)
        if not os.path.exists(TFhtml):
            os.mkdir(TFhtml)

        # one dict for all main html files. motifs_*_html_dict will contain a list of motif_page_html_dicts
        html_dict = {}
        html_dict['organism'] = organism

        bn = os.path.join(outdir, TFname)

        # 1:
        Z_hist = bn + '_FgBg-peakcall_newnoise/Z_hist.pdf'
        Z_revcum = bn + '_FgBg-peakmerge/revcum.png'

        regions = bn + '_FgBg-peakmerge/allpeaks'
        peakmerge_log = bn + '_FgBg-peakmerge/PeakMerger_log'
        peaks = bn + '_FgBg-recompz/allpeaks'
        idfile = bn + '_FgBg-peakmerge/IDfile'
        plotdir = bn + '_FgBg-selectpeaks/peak_plots.tar.gz'
        peak_plots_list = findPeakPlots(peaks, idfile, plotdir)

        os.system('convert -background white -flatten -density 90 \'%s\' \'%s/Z_hist.png\'' %(Z_hist, TFhtml)) #try to make this plot as big as matplotlib plots.
        os.system('cp \'%s\' \'%s\'' %(Z_revcum, TFhtml))

        html_dict['zhist'] = 'Z_hist.png'
        html_dict['revcum'] = 'revcum.png'

        pdir = os.path.join(TFhtml, 'peak_plots')
        os.system('mkdir \'%s\'' %pdir)


        open_archive = tarfile.open(plotdir, 'r:gz')
        for i, pl in enumerate(peak_plots_list):
            pname = os.path.split(pl)[1]
            outplot = '%i_%s' %(i+1, pname)
            open_archive.extract('./%s' %pname)
            os.system('mv \'%s\' \'%s/%s\'' %(pname, pdir, outplot))
            html_dict['peak_plot_%i' %(i+1)] = os.path.join('peak_plots', outplot)


        # for i, pl in enumerate(peak_plots_list):
        #     outplot = '%i_%s' %(i+1, os.path.split(pl)[1])
        #     os.system('cp \'%s\' \'%s/%s\'' %(pl, pdir, outplot))
        #     html_dict['peak_plot_%i' %(i+1)] = os.path.join('peak_plots', outplot)


        # PeakMerger_log
        # Z value cut-off: 4.26
        # Number of windows above Z value cut-off: 3277
        # Number of regions above Z value cut-off: 1793

        pml = open(peakmerge_log)
        z_co = pml.readline().strip().split()[-1]
        win_num = pml.readline().strip().split()[-1]
        region_num = pml.readline().strip().split()[-1]
        html_dict['z_cutoff'] = z_co
        html_dict['windows_num'] = win_num
        html_dict['regions_num'] = region_num

        peakcall_log = os.path.join(TFhtml, 'peakcall_log')
        os.system('cat \'%s\' | wc -l >> \'%s\'' %(peaks, peakcall_log))        

        nums = open(peakcall_log).readlines()
        html_dict['peaks_num'] = nums[0]
        os.system('rm \'%s\'' %peakcall_log)


        # 2:
        # count number of denovo and known WMs
        # e.g.: USF2_denovoWMstats_3-latexfrag1
        print 'Collect complementary WMs'

        d = os.listdir('%s' %outdir)

        motifs = [i for i in d if re.search(TFname+'_WMstats', i)]

        motif_nums = []
        for i in motifs:
            t = i.split('_')[-1]
            n = int(t.split('-')[0])

            motif_nums.append(n)

        motif_nums = unique(motif_nums)

        ## de novo motifs
        html_dict['motif_list'] = []

        # extract motif ensemble enrichment score (cumulative!) from enrichmentScores log:
        ensemble_enrichment_d = {}
        for line in open(bn + '_FgBg-enrichmentScores_combined_motifs/EnrichmentScores'): 
            if line.startswith('WM_path'):
                continue
            t = line.strip().split()
            name = os.path.split(t[0])[1]
            ensemble_enrichment = float(t[1])
            ensemble_enrichment_d[name] = round(ensemble_enrichment,3)

        # extract motif ensemble LL-ratio (cumulative!) from enrichmentScores log:
        ensemble_LLratio_d = {}
        for line in open(bn + '_FgBg-enrichmentScores_combined_motifs/EnrichmentScores'): 
            if line.startswith('WM_path'):
                continue
            t = line.strip().split()
            name = os.path.split(t[0])[1]
            ensemble_LLratio = float(t[2])
            ensemble_LLratio_d[name] = round(ensemble_LLratio,3)

        # read enrichment scores of separate motifs into a dictionary
        enrichment_scores_d = {}
        for line in open(bn + '_FgBg-enrichmentScores_for_all/EnrichmentScores'): 
            if line.startswith('WM_path'):
                continue
            t = line.strip().split()
            name = os.path.split(t[0])[1]
            enrichment_score = float(t[1])
            enrichment_scores_d[name] = round(enrichment_score,3)

        # read enrichment scores of separate motifs into a dictionary
        LLratio_d = {}
        for line in open(bn + '_FgBg-enrichmentScores_for_all/EnrichmentScores'): 
            if line.startswith('WM_path'):
                continue
            t = line.strip().split()
            name = os.path.split(t[0])[1]
            LLratio = float(t[3])
            LLratio_d[name] = round(LLratio,3)

        # matrix with wm as rows and columns: score, auc, fov, eabs
        motifs_dir = 'motifs'
        os.mkdir('%s' %os.path.join(TFhtml, motifs_dir))

        for wmi in motif_nums:
            wmi_html_dict = {}

            createlog = bn + '_WMstats_' + str(wmi) + '-createlogo/log_file'
            statslog = bn + '_WMstats_' + str(wmi) + '-getstats/log_file'
            identifylog = bn + '_WMstats_' + str(wmi) + '-identifymotif/top_motifs'

            f = open(createlog)
            fl = f.readlines()
            f.close()

            wmname = os.path.split(fl[0].strip().split(':')[1].lstrip())[1] #read name from original path
            wmout = os.path.join(motifs_dir, wmname)
            wmi_html_dict['name'] = urllib.quote(wmname)
            wmi_html_dict['name_unquote'] = wmname
            wmi_html_dict['score'] = enrichment_scores_d[wmname]
            wmi_html_dict['auc'] = fl[3].strip().split()[-1]
            wmi_html_dict['html'] = urllib.quote('%s.html' %wmname) #relative path from the motifs_known.html file.
            wmi_html_dict['ensemble_score'] = ensemble_enrichment_d[wmname]
            wmi_html_dict['LLratio'] = LLratio_d[wmname]
            wmi_html_dict['ensemble_LLratio'] = ensemble_LLratio_d[wmname]

            print TFname, wmout
            WMhtml = '%s/%s' %(TFhtml, wmout)
            ######
            try:
                os.system('mkdir \'%s\'' %WMhtml)
            except Exception:
                wmhtml1 = WMhtml + '_1'
                os.system('mkdir \'%s\'' %wmhtml1)
            ######

            # Copy WM that I find through _command of createlogo component
            create_command = bn + '_WMstats_' + str(wmi) + '-createlogo/_command'
            for line in open(create_command):
                if line.startswith('input.WM'):
                    wmpath = line.strip().split('=')[1]
                    new_wmpath = os.path.join(WMhtml, 'WM.txt') #wmname)
                    os.system('cp \'%s\' \'%s\'' %(wmpath, new_wmpath))
                    break

            wmi_html_dict['file'] = urllib.quote(os.path.join(wmout, 'WM.txt'))

            f = open(statslog)
            fl = f.readlines()
            f.close()

            wmi_html_dict['truepeaks'] = fl[1].strip().split()[-1]
            wmi_html_dict['eabs'] = fl[6].strip().split()[-1]
            wmi_html_dict['fov'] = fl[7].strip().split()[-1]


            wmi_html_dict['similar'] = []
            for i, fl in enumerate(open(identifylog)):
                t = fl.strip().split()
                wmi_html_dict['similar'].append({})
                wmi_html_dict['similar'][i]['name'] = t[0]
                wmi_html_dict['similar'][i]['dist'] = t[1]
                wmi_html_dict['similar'][i]['logo'] = urllib.quote('%s/%s.png' %(wmout, t[0]))


                os.system('convert \'%s/%s.pdf\' -background white -flatten \'%s/%s.png\'' %(bn + '_WMstats_' + str(wmi) + '-identifymotif/logos', t[0], WMhtml, t[0]))


            os.system('convert \'%s\' -background white -flatten \'%s/Logo.png\'' %(bn + '_WMstats_' + str(wmi) + '-createlogo/Logo.pdf', WMhtml))
            os.system('convert \'%s\' -background white -flatten \'%s/Logo_rev.png\'' %(bn + '_WMstats_' + str(wmi) + '-createlogo/Logo_rev.pdf', WMhtml))
            os.system('cp \'%s\' \'%s\'' %(bn + '_WMstats_' + str(wmi) + '-createlogo/sens_spec.png', WMhtml))

            os.system('cp \'%s\' \'%s\'' %(bn + '_WMstats_' + str(wmi) + '-getstats/zscore_post_violin.png', WMhtml))
            os.system('cp \'%s\' \'%s\'' %(bn + '_WMstats_' + str(wmi) + '-getstats/TFBS_peakcenter_dist_hist.png', WMhtml))
            os.system('cp \'%s\' \'%s\'' %(bn + '_WMstats_' + str(wmi) + '-getstats/coverage_histograms.png', WMhtml))

            os.system('cp \'%s\' \'%s\'' %(bn + '_WMstats_' + str(wmi) + '-combineWM/pdfreport.pdf', WMhtml))

            wmi_html_dict['logo'] = urllib.quote('%s/Logo.png' %wmout)
            wmi_html_dict['logo_rev'] = urllib.quote('%s/Logo_rev.png' %wmout)
            wmi_html_dict['sens_spec'] = urllib.quote('%s/sens_spec.png' %wmout)
            wmi_html_dict['zscore_post'] = urllib.quote('%s/zscore_post_violin.png' %wmout)
            wmi_html_dict['eabs_hist'] = urllib.quote('%s/coverage_histograms.png' %wmout)
            wmi_html_dict['report_pdf'] = urllib.quote('%s/pdfreport.pdf' %wmout)

            html_dict['motif_list'].append(wmi_html_dict)


            # render motif page
            motif_page_t = temp_env.get_template('motif_page.html')
            rendered_motif_page_template = motif_page_t.render(wmi_html_dict)
            o = open('%s/%s' %(TFhtml, wmname + ('.html')), 'w')
            for l in rendered_motif_page_template:
                o.write(l)
            o.close()



        # 3:
        os.system('cp \'%s\' \'%s\'' %(bn + '_motif_correlation/correlation_plot.png',os.path.join(TFhtml, motifs_dir)))

        os.system('cp \'%s\' \'%s\'' %(bn + '_FgBg-enrichmentScores_combined_motifs/EnrichmentScores.png', os.path.join(TFhtml, motifs_dir)))

        # for html:
        html_dict['motif_contribution_plot'] = os.path.join(motifs_dir, 'EnrichmentScores.png')
        html_dict['motif_correlation_heatmap'] = os.path.join(motifs_dir, 'correlation_plot.png')

        ##########################################################################
        ## Collect tops WMs stuff...
        print 'Collect tops WMs'

        d = os.listdir('%s' %outdir)

        motifs_tops = [i for i in d if re.search(TFname+'_WM_tops_stats', i)]

        motif_tops_nums = []
        for i in motifs_tops:
            t = i.split('_')[-1]
            n = int(t.split('-')[0])

            motif_tops_nums.append(n)

        motif_tops_nums = unique(motif_tops_nums)

        html_dict['motif_tops_list'] = []

        # This was alread done above:
        # motifs_dir = 'motifs'
        # os.mkdir('%s' %os.path.join(TFhtml, motifs_dir))

        for wmi in motif_tops_nums:
            wmi_html_dict = {}

            createlog = bn + '_WM_tops_stats_' + str(wmi) + '-createlogo/log_file'
            statslog = bn + '_WM_tops_stats_' + str(wmi) + '-getstats/log_file'
            identifylog = bn + '_WM_tops_stats_' + str(wmi) + '-identifymotif/top_motifs'

            f = open(createlog)
            fl = f.readlines()
            f.close()

            wmname = os.path.split(fl[0].strip().split(':')[1].lstrip())[1] #read name from original path
            wmout = os.path.join(motifs_dir, wmname)
            wmi_html_dict['name'] = urllib.quote(wmname)
            wmi_html_dict['name_unquote'] = wmname
            wmi_html_dict['score'] = round(float(fl[1].strip().split()[-1]), 3)
            wmi_html_dict['LLratio'] = round(float(fl[2].strip().split()[-1]), 3)
            wmi_html_dict['auc'] = fl[3].strip().split()[-1]
            wmi_html_dict['html'] = urllib.quote('%s.html' %wmname) #relative path from the motifs_known.html file.

            print TFname, wmout
            WMhtml = '%s/%s' %(TFhtml, wmout)
            if not os.path.exists(WMhtml):
                print 'created newly'
                os.system('mkdir \'%s\'' %WMhtml)

            # Copy WM that I find through _command of createlogo component
            create_command = bn + '_WM_tops_stats_' + str(wmi) + '-createlogo/_command'
            for line in open(create_command):
                if line.startswith('input.WM'):
                    wmpath = line.strip().split('=')[1]
                    new_wmpath = os.path.join(WMhtml, 'WM.txt') #wmname)
                    os.system('cp \'%s\' \'%s\'' %(wmpath, new_wmpath))
                    break

            wmi_html_dict['file'] = urllib.quote(os.path.join(wmout, 'WM.txt'))

            f = open(statslog)
            fl = f.readlines()
            f.close()

            wmi_html_dict['truepeaks'] = fl[1].strip().split()[-1]
            wmi_html_dict['eabs'] = fl[6].strip().split()[-1]
            wmi_html_dict['fov'] = fl[7].strip().split()[-1]


            wmi_html_dict['similar'] = []
            for i, fl in enumerate(open(identifylog)):
                t = fl.strip().split()
                wmi_html_dict['similar'].append({})
                wmi_html_dict['similar'][i]['name'] = t[0]
                wmi_html_dict['similar'][i]['dist'] = t[1]
                wmi_html_dict['similar'][i]['logo'] = urllib.quote('%s/%s.png' %(wmout, t[0]))


                os.system('convert \'%s/%s.pdf\' -background white -flatten \'%s/%s.png\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-identifymotif/logos', t[0], WMhtml, t[0]))


            os.system('convert \'%s\' -background white -flatten \'%s/Logo.png\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-createlogo/Logo.pdf', WMhtml))
            os.system('convert \'%s\' -background white -flatten \'%s/Logo_rev.png\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-createlogo/Logo_rev.pdf', WMhtml))
            os.system('cp \'%s\' \'%s\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-createlogo/sens_spec.png', WMhtml))

            os.system('cp \'%s\' \'%s\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-getstats/zscore_post_violin.png', WMhtml))
            os.system('cp \'%s\' \'%s\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-getstats/TFBS_peakcenter_dist_hist.png', WMhtml))
            os.system('cp \'%s\' \'%s\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-getstats/coverage_histograms.png', WMhtml))

            os.system('cp \'%s\' \'%s\'' %(bn + '_WM_tops_stats_' + str(wmi) + '-combineWM/pdfreport.pdf', WMhtml))

            wmi_html_dict['logo'] = urllib.quote('%s/Logo.png' %wmout)
            wmi_html_dict['logo_rev'] = urllib.quote('%s/Logo_rev.png' %wmout)
            wmi_html_dict['sens_spec'] = urllib.quote('%s/sens_spec.png' %wmout)
            wmi_html_dict['zscore_post'] = urllib.quote('%s/zscore_post_violin.png' %wmout)
            wmi_html_dict['eabs_hist'] = urllib.quote('%s/coverage_histograms.png' %wmout)
            wmi_html_dict['report_pdf'] = urllib.quote('%s/pdfreport.pdf' %wmout)

            html_dict['motif_tops_list'].append(wmi_html_dict)


            # render motif page
            motif_page_t = temp_env.get_template('motif_page.html')
            rendered_motif_page_template = motif_page_t.render(wmi_html_dict)
            o = open('%s/%s' %(TFhtml, wmname + ('.html')), 'w')
            for l in rendered_motif_page_template:
                o.write(l)
            o.close()


        ##########################################################################


        # 4:
        os.system('cp \'%s\' \'%s\'' %(bn + '_FgBg-combine2/pdfreport.pdf', TFhtml))
        html_dict['report_pdf'] = 'pdfreport.pdf'

        # 5:
        os.system('cp \'%s\' \'%s/peaks_with_sites.txt\'' %(bn + '_combinepeakstfbs/peaks_with_sites', TFhtml))

        html_dict['peak_lines'] = []
        motif_names = []
        for i, line in enumerate(open(os.path.join(TFhtml, 'peaks_with_sites.txt'))):
            if i == 0:
                #get all motif names:
                if organism != 'dm3':
                    t = line.split('(##-separated)')[1]
                else:
                    t = line.split('gene_strand')[-1]

                t = t.split('binding_sites')[0]
                t = t.split()
                for n in t:
                    motif_names.append(re.sub('_site_number$', '', n))
                html_dict['motif_names'] = motif_names
                continue

            line_dict = {}
            t = line.strip().split('\t')
            line_dict['chrom'] = t[0]
            line_dict['start'] = int(t[1])
            line_dict['stop'] = int(t[2])
            line_dict['score'] = t[4]
            # t[6] are nearest upstream promoter members and t[7] its offset. t[12] and t[13] the same for downstream. t[18] are gene descriptions.
            if organism != 'dm3':
                line_dict['upstream_genes'] = parseGenes(t[6], t[24])
                line_dict['upstream_offset'] = t[7] + ' (%s)' %t[8]
                line_dict['downstream_genes'] = parseGenes(t[15], t[24])
                line_dict['downstream_offset'] = t[16] + ' (%s)' %t[17]
            else:
                #chrom  start   stop    peakID  Z-score strand  nearest_upstream_gene_name|ensembl_id   offset  gene_strand     2nearest_upstream_gene_name|ensembl_id  offset  gene_strand     3nearest_upstream_gene_name|ensembl_id  offset  gene_strand     nearest_downstream_gene_name|ensembl_id offset  gene_strand     2nearest_downstream_gene_name|ensembl_id        offset  gene_strand     3nearest_downstream_gene_name|ensembl_id        offset  gene_strand    denovo_WM_12_site_number        denovo_WM_10_site_number        Cg12029_site_number     Z_site_number   binding_sites (WM_name:relative distance to peak center:posterior binding probability)
                line_dict['upstream_genes'] = t[6]
                line_dict['upstream_offset'] = t[7] + ' (%s)' %t[8]
                line_dict['downstream_genes'] = t[15]
                line_dict['downstream_offset'] = t[16] + ' (%s)' %t[17]

            # adapt the frame of the regions shown in swiss regulon browser
            try:
                if abs(int(t[7])) < 10000 or abs(int(t[13])) < 10000:
                    if abs(int(t[7])) <= abs(int(t[13])):
                        swiss_reg_frame_start = int(t[1]) + int(t[7]) - 1000
                        swiss_reg_frame_stop = int(t[2]) + 1000
                    else:
                        swiss_reg_frame_start = int(t[1]) - 1000
                        swiss_reg_frame_stop = int(t[2]) + int(t[13]) + 1000
                else:
                    swiss_reg_frame_start = int(t[1]) - 1000
                    swiss_reg_frame_stop = int(t[2]) + 1000
            except ValueError:
                if t[7] == 'None' and t[13] == 'None':
                    swiss_reg_frame_start = int(t[1]) - 1000
                    swiss_reg_frame_stop = int(t[2]) + 1000
                else:
                    if t[7] != 'None':
                        if abs(int(t[7])) < 10000:
                            swiss_reg_frame_start = int(t[1]) + int(t[7]) - 1000
                            swiss_reg_frame_stop = int(t[2]) + 1000
                    else:
                        if abs(int(t[13])) < 10000:
                            swiss_reg_frame_start = int(t[1]) - 1000
                            swiss_reg_frame_stop = int(t[2]) + int(t[13]) + 1000

            line_dict['swiss_reg_frame_start'] = swiss_reg_frame_start
            line_dict['swiss_reg_frame_stop'] = swiss_reg_frame_stop


            # binding sites are at the end after gene descriptions t[18] and site_numbers 
            line_dict['sites'] = getSitesDict(motif_names, t[24+len(motif_names)+1:])
            html_dict['peak_lines'].append(line_dict)
            if i == 200:
                break

        html_dict['peak_file'] = 'peaks_with_sites.txt'


        # 6: add replicate stuff (wig, bed, pdf) to download
        html_dict['samples'] = TFs_samples_html_dict[TFname]

        # add background sample dicts
        for d in TFs_samples_html_dict['BG']:
            html_dict['samples'].append(d)

        # add project name
        html_dict['project'] = ''
        logging.warning(os.getcwd())
        if os.path.exists('project'):
            logging.warning("adding project")
            html_dict['project'] = open("project").read().split("\n")[0].rstrip()

        # render html
        html_temp = temp_env.get_template('index.html')
        rendered_html_temp = html_temp.render(html_dict)
        o = open(os.path.join(TFhtml,'index.html'), 'w')
        for l in rendered_html_temp:
            o.write(l)
        o.close()


        # copy css and js directories
        os.system('cp -r %s %s' %(os.path.join(pipeline_dir, 'templates/css') , TFhtml))
        os.system('cp -r %s %s' %(os.path.join(pipeline_dir, 'templates/js') , TFhtml))
        os.system('cp -r %s %s' %(os.path.join(pipeline_dir, 'templates/images') , TFhtml))


if __name__ == '__main__':
    main()
