#!/usr/bin/awk -f

BEGIN {
    # Parameters that you will set via -v flags on the command line
    # For example: -v GAP_LENGTH=10 -v MIN_LENGTH=100
    GAP_LENGTH   = 10;     # Maximum gap between positions in the same region
    MIN_LENGTH   = 100;    # Minimum length for a region to be reported

    # G1_NUM and G2_NUM will also be passed via -v, e.g., -v G1_NUM=5 -v G2_NUM=7
    # G1_NUM: Number of samples in Group 1
    # G2_NUM: Number of samples in Group 2
}

{
    # Each line represents a genomic position with depth data for all samples.
    # $1: Contig ID
    # $2: Genomic Position
    # $3 onwards: Depth values for samples

    ctg = $1;
    pos = $2;

    # Flags to check if ALL samples in G1 and ALL samples in G2 have depth > 0 at this position
    all_G1_depth_gt_0 = 1; # Assume true initially
    all_G2_depth_gt_0 = 1; # Assume true initially

    # --- Check ALL G1 samples ---
    # G1 samples start at column 3 ($3)
    # They go up to column (2 + G1_NUM)
    # Example: If G1_NUM=3, G1 samples are $3, $4, $5. The loop goes from i=3 to (2+3)=5.
    for (i = 3; i <= 2 + G1_NUM; i++) {
        if ($i <= 0) {
            all_G1_depth_gt_0 = 0; # Found a G1 sample with non-positive depth
            break;                 # No need to check other G1 samples for this position
        }
    }

    # --- Check ALL G2 samples ---
    # G2 samples start immediately after the last G1 sample.
    # The first G2 sample is at column (3 + G1_NUM)
    # They go up to the last column (NF)
    # Example: If G1_NUM=3, G2 samples start at (3+3)=6. They go from $6 to $NF.
    for (i = 3 + G1_NUM; i <= NF; i++) {
        if ($i <= 0) {
            all_G2_depth_gt_0 = 0; # Found a G2 sample with non-positive depth
            break;                 # No need to check other G2 samples for this position
        }
    }

    # --- Combine conditions: Position is valid only if ALL G1 AND ALL G2 samples have depth > 0 ---
    if (all_G1_depth_gt_0 && all_G2_depth_gt_0) {
        # Store this position for later grouping into regions
        # We use a space as a separator to easily split the positions later
        BothGT0_CTG[ctg] = BothGT0_CTG[ctg] " " pos;
    }
}

END {
    # This block runs after all lines in the input file have been processed.

    # --- Group the identified positions into continuous regions ---
    # Iterate through each contig that had at least one valid position
    for (ctg in BothGT0_CTG) {
        # Split the space-separated positions string into an array
        split(BothGT0_CTG[ctg], pos_arr, " ");

        # Sort the positions numerically. This is crucial for correctly identifying regions.
        # 'asort' is a GNU AWK (gawk) extension for numerical sorting.
        # If you're using a very basic AWK, this might not work, and you'd need
        # to ensure your input positions are already sorted or use an external sort.
        asort(pos_arr);

        start = 0; # Start position of the current potential region
        last = 0;  # Last position included in the current potential region
        end = 0;   # Current position being evaluated

        # Loop through the sorted positions for the current contig
        for (i = 1; i <= length(pos_arr); i++) {
            end = pos_arr[i]; # Get the current position

            if (i == 1) {
                # This is the very first position for this contig, so it starts a new region.
                start = end;
                last = end;
                continue; # Move to the next position
            }

            # Check if the current position is within the allowed gap from the last position
            if (end - last <= GAP_LENGTH) {
                last = end; # Extend the current region
            } else {
                # Gap is too large, the current region has ended.
                # Calculate its length (inclusive: last - start + 1)
                len = last - start + 1;

                # If the region meets the minimum length requirement, print it
                if (len >= MIN_LENGTH) {
                    # Output to a file for regions where both G1 and G2 have depth > 0
                    print ctg "\t" start "\t" last "\t" len > "G1_and_G2_gt0_regions.txt";
                }
                
                # Start a new region with the current 'end' position
                start = end;
                last = end;
            }
        }
        # After the loop finishes, there might be a pending region.
        # This part ensures the very last identified region is processed.
        len = last - start + 1;
        if (len >= MIN_LENGTH) {
            print ctg "\t" start "\t" last "\t" len > "G1_and_G2_gt0_regions.txt";
        }
    }
    # Close the output file after all processing is done
    close("G1_and_G2_gt0_regions.txt");
}