awk -F'\t' 'NR > 1 {print $7}' /home/phuctran/Desktop/1.clean_eggnog.annotations | sort | uniq -c > COG.tsv

awk '{
    # Skip the header line (adjust based on your file)
    if ($2 == "-") { next }

    # Get the count and the category string
    count = $1;
    category = $2;

    # Loop through each character (letter) in the category string
    for (i = 1; i <= length(category); i++) {
        
        # Extract the single letter
        letter = substr(category, i, 1);
        
        # Add the count to the total for that single letter
        counts[letter] += count;
    }
}
END {
    # Print the final tallied counts
    for (letter in counts) {
        print letter, counts[letter];
    }
}' COG.tsv | sort > COG_count.tsv

