#include <stdio.h>

// Function to print an array
void printArray(int arr[], int n) {
    for (int i = 0; i < n; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

void heap(int arr[], int size, int n) {
    if (size == 1) {
        printArray(arr, n);
        return;
    }

    for (int i = 0; i < size; i++) {
        heap(arr, size - 1, n);

        if (size % 2 == 1) {
            int temp = arr[0];
            arr[0] = arr[size - 1];
            arr[size - 1] = temp;
        } else {
            int temp = arr[i];
            arr[i] = arr[size - 1];
            arr[size - 1] = temp;
        }
    }
}

void 
// Driver code
int main() {
    int arr[] = {1, 2, 3, 4};

    printf("All permutations:\n");
    heap(arr, 3, 4);

    return 0;
}