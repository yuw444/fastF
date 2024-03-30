#include "data.h"

umi_node *new_umi_node(char *umi)
{
    if (umi == NULL) {
        return NULL;
    }

    umi_node *node = (umi_node *)malloc(sizeof(umi_node));
    node->umi = (char *)malloc(strlen(umi) + 1);
    strcpy(node->umi, umi);
    node->left = NULL;
    node->right = NULL;
    return node;
}


umi_node *insert_umi_node(umi_node *root, char *umi)
{
    if (root == NULL) {
        return new_umi_node(umi);
    }

    if (strcmp(umi, root->umi) < 0) {
        root->left = insert_umi_node(root->left, umi);
    } else if (strcmp(umi, root->umi) > 0) {
        root->right = insert_umi_node(root->right, umi);
    }

    return root;
}

void free_umi_node(umi_node *root)
{
    if (root == NULL) {
        return;
    }

    free_umi_node(root->left);
    free_umi_node(root->right);
    free(root->umi);
    free(root);
}


gene_umi_node *new_gene_umi_node(char *gene, char *umi)
{
    if (gene == NULL || umi == NULL) {
        return NULL;
    }

    gene_umi_node *node = (gene_umi_node *)malloc(sizeof(gene_umi_node));
    node->gene = (char *)malloc(strlen(gene) + 1);
    strcpy(node->gene, gene);
    node->nth_gene = 0;
    node->expression_level = 0;
    node->umi_root = new_umi_node(umi);
    node->left = NULL;
    node->right = NULL;
    return node;
}

gene_umi_node *insert_gene_umi_node(gene_umi_node *root, char *gene, char *umi)
{
    if (root == NULL) {
        return new_gene_umi_node(gene, umi);
    }

    if (strcmp(gene, root->gene) < 0) {
        root->left = insert_gene_umi_node(root->left, gene, umi);
    } else if (strcmp(gene, root->gene) > 0) {
        root->right = insert_gene_umi_node(root->right, gene, umi);
    } else {
        root->umi_root = insert_umi_node(root->umi_root, umi);
    }

    return root;
}

void free_gene_umi_node(gene_umi_node *root)
{
    if (root == NULL) {
        return;
    }

    free_gene_umi_node(root->left);
    free_gene_umi_node(root->right);
    free(root->gene);
    free_umi_node(root->umi_root);
    free(root);
}

