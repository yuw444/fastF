
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct data_node {
    int data;
    struct data_node* left;
    struct data_node* right;
}data_node;

data_node* new_data_node(int data) {
    data_node* node = (data_node*)malloc(sizeof(data_node));
    node->data = data;
    node->left = NULL;
    node->right = NULL;
    return node;
}

data_node* insert_data_node(data_node* root, int data) {
    if (root == NULL) {
        return new_data_node(data);
    }
    
    if (data < root->data) {
        root->left = insert_data_node(root->left, data);
    }
    else if (data > root->data) {
        root->right = insert_data_node(root->right, data);
    }
    
    return root;
}

void inorder(data_node* root) {
    if (root == NULL) {
        return;
    }
    
    inorder(root->left);
    printf("%d ", root->data);
    inorder(root->right);
}


void inorder_index(data_node* root, int data, int* index) {
    if (root == NULL) {
        return;
    }
    
    inorder_index(root->left, data, index);
    if (root->data <= data) {
        (*index)++;
    }
    inorder_index(root->right, data, index);
}


int main(int argc, char *argv[])
{

    data_node *root = new_data_node(5);
    root = insert_data_node(root, 3);
    root = insert_data_node(root, 4);
    root = insert_data_node(root, 7);
    root = insert_data_node(root, 6);
    root = insert_data_node(root, 8);
    root = insert_data_node(root, 1);

    inorder(root);

    int index = 0;
    for (int i = 1; i < 9; i++) {
        inorder_index(root, i, &index);
        printf("order of : %d is %d \n", i, index);
        index = 0;
    }
    return 0;

}
