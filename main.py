import asyncio
import aiohttp
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Data import CodonTable
import os
import hashlib
import io
import re
import asyncio


file_path = "example.fasta"
output_file = "blast_results.txt"

# Get the standard genetic code table
genetic_code = CodonTable.unambiguous_dna_by_name["Standard"]

# Initialize the BLAST results cache
blast_results_cache = {}


async def blast_search(protein_sequence):
    # Check if the result is already in the cache
    sequence_hash = hashlib.md5(str(protein_sequence).encode()).hexdigest()
    if sequence_hash in blast_results_cache:
        print('Got one result from cache')
        return blast_results_cache[sequence_hash]

    # Perform online BLAST search
    async with aiohttp.ClientSession() as session:
        # Step 1: Submit the BLAST query and get the request ID
        async with session.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data={
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "swissprot",
            "QUERY": str(protein_sequence)
        }) as response:
            response_text = await response.text()
            request_id = re.search(r"RID = (\w+)", response_text).group(1)

        # Step 2: Check the status of the BLAST query and wait for it to complete
        status = "WAITING"
        while status in ["WAITING", "QUEUED"]:
            await asyncio.sleep(10)  # Wait for 10 seconds before checking the status again
            async with session.get(f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={request_id}") as status_response:
                status_text = await status_response.text()
                status = re.search(r"Status=(\w+)", status_text).group(1)

        # Step 3: Once the query is done, fetch the BLAST results in XML format
        if status == "READY":
            print('BLAST query is ready')
            async with session.get(f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={request_id}") as result_response:
                result_text = await result_response.text()
                if result_text.startswith('<?xml'):
                    blast_record = NCBIXML.read(io.StringIO(result_text))
                    # Cache the result
                    blast_results_cache[sequence_hash] = blast_record
                    return blast_results_cache[sequence_hash]
                else:
                    raise ValueError("BLAST response was not in XML format.")
        else:
            raise ValueError(f"BLAST query failed with status: {status}")


async def process_sequences(sequences):
    tasks = []
    for seq_record in sequences:
        dna_sequence = seq_record.seq
        # Trim sequence to make it a multiple of three
        trimmed_sequence = dna_sequence[:len(dna_sequence) - (len(dna_sequence) % 3)]
        # Translate all six frames and create BLAST tasks
        for frame in range(6):
            frame_sequence = trimmed_sequence[frame:]
            protein_sequence = frame_sequence.translate(table=genetic_code, to_stop=False)
            tasks.append(blast_search(protein_sequence))
            print('One Frame done')
    return await asyncio.gather(*tasks)
try:
    with open(file_path, "r") as file, open(output_file, "w") as output:
        sequences = SeqIO.parse(file, "fasta")
        results = asyncio.run(process_sequences(sequences))
        # Write results to the output file
        for blast_record in results:
            for alignment in blast_record.alignments:
                output.write(f"\nMatch: {alignment.title}\n")
                for hsp in alignment.hsps:
                    output.write(f"Bit score: {hsp.bits}\n")
                    output.write(f"E-value: {hsp.expect}\n")
                    output.write(f"Alignment length: {hsp.align_length}\n")
                    output.write(f"Query sequence: {hsp.query}\n")
                    output.write(f"Match sequence: {hsp.match}\n")
                    output.write(f"Subject sequence: {hsp.sbjct}\n")
    print(f"BLAST results saved to {output_file}")
except IOError as e:
    print(e)
    # Create the output file if it does not exist
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_file, "w") as output:
        output.write("BLAST results will be saved here.")