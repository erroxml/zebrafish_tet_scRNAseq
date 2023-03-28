cellranger count --id=WT1 \
                   --transcriptome=/home/wll/genomes/GRCz11/GRCz11_cellranger/Danio.rerio_genome \
                   --fastqs=/home/wll/projects/scTet20220707/01raw_data/WT1 \
                   --sample=WT1 \
                   --expect-cells=10000 \
                   --localcores=24 \
                   --localmem=100
cellranger count --id=WT2 \
                   --transcriptome=/home/wll/genomes/GRCz11/GRCz11_cellranger/Danio.rerio_genome \
                   --fastqs=/home/wll/projects/scTet20220707/01raw_data/WT2 \
                   --sample=WT2 \
                   --expect-cells=10000 \
                   --localcores=24 \
                   --localmem=100
cellranger count --id=KO1 \
                   --transcriptome=/home/wll/genomes/GRCz11/GRCz11_cellranger/Danio.rerio_genome \
                   --fastqs=/home/wll/projects/scTet20220707/01raw_data/KO1 \
                   --sample=KO1 \
                   --expect-cells=10000 \
                   --localcores=24 \
                   --localmem=100
cellranger count --id=KO2 \
                   --transcriptome=/home/wll/genomes/GRCz11/GRCz11_cellranger/Danio.rerio_genome \
                   --fastqs=/home/wll/projects/scTet20220707/01raw_data/KO2 \
                   --sample=KO2 \
                   --expect-cells=10000 \
                   --localcores=24 \
                   --localmem=100