#include <stdio.h>
#include <stdlib.h>

// Define the structure for a binary search tree node
struct TreeNode {
    int key;
    int value;
    struct TreeNode* left;
    struct TreeNode* right;
};

// Define the structure for a hash table node
struct HashNode {
    int key;
    int value;
    struct HashNode* next;
};

// Define the hash table size (can be adjusted based on requirements)
#define HASH_TABLE_SIZE 100

// Helper function to create a new binary search tree node
struct TreeNode* newTreeNode(int key, int value) {
    struct TreeNode* node = (struct TreeNode*)malloc(sizeof(struct TreeNode));
    node->key = key;
    node->value = value;
    node->left = NULL;
    node->right = NULL;
    return node;
}

// Helper function to insert a key-value pair into the hash table
void insertHashNode(struct HashNode** hashTable, int key, int value) {
    // Compute the hash value
    int hashValue = key % HASH_TABLE_SIZE;
    // Create a new hash node
    struct HashNode* newNode = (struct HashNode*)malloc(sizeof(struct HashNode));
    newNode->key = key;
    newNode->value = value;
    newNode->next = NULL;
    // Insert the new node at the head of the linked list for the hash value
    newNode->next = hashTable[hashValue];
    hashTable[hashValue] = newNode;
}

// Helper function to convert a binary search tree to a hash table
void convertBSTtoHashTable(struct TreeNode* root, struct HashNode** hashTable) {
    if (root == NULL) {
        return;
    }
    // Traverse the left subtree
    convertBSTtoHashTable(root->left, hashTable);
    // Insert the current node's key-value pair into the hash table
    insertHashNode(hashTable, root->key, root->value);
    // Traverse the right subtree
    convertBSTtoHashTable(root->right, hashTable);
}

// free the tree
void free_tree(struct TreeNode *root)
{
    if (root == NULL)
    {
        return;
    }
    free_tree(root->left);
    free_tree(root->right);
    free(root);
}

// free the hash table
void free_hash_table(struct HashNode **hashTable)
{
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        struct HashNode* currentNode = hashTable[i];
        while (currentNode != NULL) {
            struct HashNode* temp = currentNode;
            currentNode = currentNode->next;
            free(temp);
        }
    }
}

int main() {
    // Create a sample binary search tree
    struct TreeNode* root = newTreeNode(4, 40);
    root->left = newTreeNode(2, 20);
    root->right = newTreeNode(6, 60);
    root->left->left = newTreeNode(1, 10);
    root->left->right = newTreeNode(3, 30);
    root->right->left = newTreeNode(5, 50);
    root->right->right = newTreeNode(7, 70);
    // Create a hash table array
    struct HashNode* hashTable[HASH_TABLE_SIZE] = {NULL};
    // Convert the BST to the hash table
    convertBSTtoHashTable(root, hashTable);
    // Print the hash table
    for (int i = 0; i < HASH_TABLE_SIZE; i++) {
        printf("Hash table entry %d:\n", i);
        struct HashNode* currentNode = hashTable[i];
        while (currentNode != NULL) {
            printf("Key: %d, Value: %d\n", currentNode->key, currentNode->value);
            currentNode = currentNode->next;
        }
    }

    // free the tree
    free_tree(root);
    free_hash_table(hashTable);

    return 0;
}
