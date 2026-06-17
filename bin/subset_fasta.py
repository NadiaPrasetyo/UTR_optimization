import argparse
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(input_file, output_file, accessions):
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        if record.id.split(" ")[0] in accessions:
            records.append(record)

    SeqIO.write(records, output_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("--accessions", nargs="+", help="Accessions to subset")
    args = parser.parse_args()

    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    if not args.accessions:
        accessions = ["MH806866.1","NC_023858.1","KF961188.1","NC_023858.1","MW684814.1","OR590795.1","KF979335.1","NC_024768.1","MH125198.1","PP228891.1","PP228891.1","PP228890.1","PP228897.1","KF961187.1","NC_024769.1","KF979336.1","NC_024769.1","MG846465.1","PX095761.1","NC_038957.1","KC811837.1","OR364988.1","KC811837.1","NC_038957.1","NC_039235.1","KF961186.1","NC_039235.1","OP169446.1","OP169446.1","NC_024769.1","MZ210013.1","PQ360974.1","OR837777.1","MZ436972.1","ON949942.1","PX392724.1","PX464736.1","PX580375.1","PX580376.1","PV134467.1","PQ279224.1","PQ279225.1","PQ493407.1","PQ493408.1","PQ623186.1","PP437939.1","PP437940.1","OR826802.1","OQ162219.1","ON949938.1","ON125123.1","OM108481.1","OM108482.1","OL630964.1","OK334537.1","OK334538.1",NC.034206, NC.034973, NC.012531, NC.012800, "MZ821772.1","MZ821781.1","MZ821786.1",NC.003003, NC.003113, "MW314665.1","MW314666.1","MT840202.1","MW019503.1","MN933901.1","MN423334.1","MN598012.1","MN598019.1","MN598020.1","MN598023.1","MN598027.1","MN598040.1","MN598041.1","MN663160.1","MH745407.1","MK606538.1","MK639928.1","MH201082.1","MH587229.1","MF797911.1","LC898765.1","LC898766.1","LC898794.1","KY888273.1","KY670597.1","KX068679.1","KJ789136.1","KJ789137.1","KM986843.1","KP036483.1","JX976770.1","JQ639383.1","JQ639384.1","JX017380.1","FJ496339.1","FJ496340.1","FJ496341.1","FJ496342.1","FJ496343.1","FJ496344.1","FJ528584.1","GQ485310.1","GQ485311.1","GU017972.1","HM245928.1","HQ891923.1","BK061329.1","DQ401688.1","EF428566.1","AF225953.1","AF323747.1","PQ736759.1","NC_038767.1","KC904083.1","NC_028964.1","KT880666.1","MT127560.1","MH998217.1","LC851702.1","MK834286.1","PQ238660.1","OQ210941.1","JX976773.1","LC209440.1","LC209462.1","LC209469.1","LC209465.1","LC209473.1","KY888272.1","OP093966.1","OP093968.1","BK059624.1"]
    else:
        accessions = args.accessions
    main(args.input, args.output, accessions)