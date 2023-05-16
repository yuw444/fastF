
SELECT cell_index, feature_index, COUNT(DISTINCT encoded_umi) AS umi_count 
FROM umi
GROUP BY cell_index, feature_index;