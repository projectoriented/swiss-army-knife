configfile: 'browser_snapshots.yaml'

bg=config['bg']
sessions_txt=config['sessions_txt']
make_bed_script=config['make_bed_script']
sessions_script=config['sessions_script']

rule subset_sessions:
    input:
        sessions_file = sessions_txt
    output:
        subsetted_sessions = '{sample}/{sample}-sv_sessions.txt'
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        sed -r "s/(^hub_[0-9]{{3}}_)({wildcards.sample})(_\S+)(_\S+)(\s)([a-z]+)/\\1\\2\\3\\4 squish/" {input.sessions_file} > {output.subsetted_sessions} 
        #sed -r "s/(^hub_[0-9]{{3}}_)(\S+)(_\S+)(_\S+)(\s)([a-z]+)/\\1\\2\\3\\4 squish/" {input.sessions_file} > {output.subsetted_sessions} 
        '''

rule subset_bg_tsv:
    input:
        bg_tsv = bg
    output:
        subset_bg = '{sample}/{sample}-bg.tsv'
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        param="{wildcards.sample} {input.bg_tsv}"
        if grep -q $param; then
            grep $param > {output.subset_bg}
        else
            touch {output.subset_bg}
        fi
        '''

rule get_bed:
    input:
        tsv = rules.subset_bg_tsv.output.subset_bg,
        sample_session = rules.subset_sessions.output.subsetted_sessions
    output:
        bed = '{sample}/{sample}-browser.bed'
    params:
        make_bed = make_bed_script
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        cat {input.tsv} | cut -f7 | while read line; do {params.make_bed} $line $(basename {input.sample_session}) ; done > {output.bed}
        '''

rule browser_snapshots:
    input:
        bed = rules.get_bed.output.bed
    output:
        pdf = '{sample}.pdf'
    params:
        sessions_py = sessions_script
    threads: 1
    resources:
        mem = 1,
        hrs=12
    shell:
        '''
        if [ -s {input.bed} ]; then
            cd $( dirname {input.bed} )
            {params.sessions_py} --regions_file $(basename {input.bed}) && pdfunite *.pdf {output.pdf} && mv {output.pdf} ../
        else
            touch {output.pdf}
        fi
        '''

onsuccess:
    shell('ls -d */ | xargs rm -rf')