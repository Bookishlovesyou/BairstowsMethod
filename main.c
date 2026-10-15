#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef enum Operator {
    ADD,
    SUBSTRACT
} Operator;

typedef struct Polynomial {
    unsigned short degree;
    float *terms;
}Polynomial;

typedef struct Table {
    Polynomial *first_func;
    Polynomial *middle_func;
    Polynomial *last_func;
    float neg_u0;
    float neg_v0;
    float delta_u;
    float delta_v;
    float u_k1;
    float v_k1;
}Table;

typedef struct Roots {
    float real;
    float imag;
}Roots;

typedef struct Sheet {
    float tolerance;
    Roots *roots[5];
    Polynomial *func;
    Table **tables;
}Sheet;



Polynomial *create_equation(unsigned short degree);
Sheet *create_sheet(unsigned short degree);
Table *create_table(float u0, float v0, Polynomial *first_func);

void display(Polynomial *func);

unsigned short input_degree() {
    unsigned short degree;
    printf("Enter degree of the Polynomial: ");
    scanf("%hd", &degree);
    if(degree < 2 || degree > 5) {
        printf("Degree must be greater than 1 and Less than 6!\n");
        return input_degree();
    }
    return ++degree;
}

int main() {
    unsigned short degree = input_degree();
    Sheet *sheet1 = create_sheet(degree);
    // free(sheet1);
    return 0;
}

int check_negative(float a, float b, float c) {
    return (pow(b,2) - (4 * a * c)) < 0;
}

float get_i_num(float a, float b, float c) {
    return sqrt(fabsf(pow(b,2) - (4 * a * c))) / (2 * a);
}

float get_real_num(float a, float b) {
    return ((-1) * b) / (2 * a);
}


void get_quadratic_root(Roots *root, Polynomial *func, Operator op) {

    float a = func->terms[0];
    float b = func->terms[1];
    float c = func->terms[2];


    if(!check_negative(a,b,c)) {
        switch(op) {
            case ADD:
                root->real = get_real_num(a,b) + get_i_num(a,b,c);
                break;
            case SUBSTRACT:
                root->real = get_real_num(a,b) - get_i_num(a,b,c);
                break;
        }
        
        root->imag = 0;
    }
    else {
        root->real = get_real_num(a,b);
        root->imag = get_i_num(a,b,c);
    }
}


Polynomial *create_quadtratic(float a, float b, float c) {
    Polynomial *equation = (Polynomial*)malloc(sizeof(Polynomial));
    equation->terms = (float*)malloc(3 * sizeof(float));
        equation->degree = 3;

    equation->terms[0] = a;
    equation->terms[1] = b;
    equation->terms[2] = c;
    return equation;
}

Polynomial *fetch_left_func(Polynomial *func) {
    Polynomial *new_func = (Polynomial*)malloc(sizeof(Polynomial));
    new_func->terms = (float*) malloc((func->degree-1) * sizeof(float));
    new_func->degree = func->degree-1;
    for(int i=0;i<=func->degree-1; i++) {
        new_func->terms[i] = func->terms[i];
    }
    return new_func;
}

Roots *shift_solve(Polynomial *func) {
    Roots *root = (Roots*)malloc(sizeof(Roots));
    root->real = ((-1) * func->terms[1]) / func->terms[0];
    root->imag = 0;
    return root;
}

void loop(Sheet *sheet, Polynomial *func, float tolerance) {
    unsigned short max = func->degree-1;
    Table *table = *(sheet->tables);


    while(fabsf(table->middle_func->terms[max]) + fabsf(table->middle_func->terms[max-1]) > tolerance) {
        table = create_table(table->u_k1, table->v_k1, func);
    }   
    *sheet->tables = table;
}

void tolerance_loop(Sheet *sheet, Polynomial *func, float tolerance) {
    for(int i=0;i <5;i++) {
        sheet->roots[i] = (Roots*)malloc(sizeof(Roots));
    }
    if(func->degree == 3) {
        get_quadratic_root(sheet->roots[0], func, ADD);
        get_quadratic_root(sheet->roots[1], func, SUBSTRACT);
        return;
    }
    loop(sheet, func, tolerance);
    Table *table = *sheet->tables;
    Polynomial *quadratic = create_quadtratic(1, table->u_k1, table->v_k1);
    
    if (quadratic == NULL) {
        // Handle error
        return;
    }
    
    Polynomial *shorten_func = fetch_left_func(table->middle_func);
    
    if (quadratic && shorten_func) {


        get_quadratic_root(sheet->roots[0],quadratic, ADD);
        get_quadratic_root(sheet->roots[1],quadratic, SUBSTRACT);
        
        switch(func->degree) {
            case 4:
                sheet->roots[2] = shift_solve(shorten_func);
                break;
            case 5:
                get_quadratic_root(sheet->roots[2], table->middle_func, ADD);
                get_quadratic_root(sheet->roots[3], table->middle_func, SUBSTRACT);
                break;
            case 6:
                shorten_func->degree = func->degree - 2;
                *sheet->tables = create_table(table->u_k1, table->v_k1, shorten_func);
                loop(sheet, shorten_func, tolerance);
                table = *sheet->tables;
                if (table) {
                    Polynomial *quadratic2 = create_quadtratic(1, table->u_k1, table->v_k1);
                    if (quadratic2) {
                        get_quadratic_root(sheet->roots[2], quadratic2, ADD);
                        get_quadratic_root(sheet->roots[3], quadratic2, SUBSTRACT);
                        sheet->roots[4] = shift_solve(fetch_left_func(table->middle_func));
                        free(quadratic2);
                    }
                }
                break;
        }
    }
    
    if (quadratic) {
        free(quadratic->terms);
        free(quadratic);
    }
    if (shorten_func) {
        free(shorten_func->terms);
        free(shorten_func);
    }
}

void display_roots(Roots **roots, unsigned short degree) {
    printf("Roots:\n");

    for(int i=0; i<degree-1; i++) {
        if(roots[i]->imag != 0) {
            printf("x%d = %.6f %c %.6fi\n", i+1, roots[i]->real,(i % 2 == 0) ? '+' : '-', roots[i]->imag);
        }
        else {
            printf("x%d = %.6f\n",i+1, roots[i]->real);
        }
    }
    printf("\n");
}


Sheet *create_sheet(unsigned short degree) {
    Sheet *paper = (Sheet*) malloc(sizeof(Sheet));
    paper->tables = (Table**)malloc(sizeof(Table*));
    paper->func = create_equation(degree);
    // paper->roots = (Roots**)malloc((degree+1) * sizeof(Roots*));
    printf("Enter Tolerance: ");
    scanf("%f", &paper->tolerance);
    *paper->tables = create_table(1,1, paper->func);
    tolerance_loop(paper, paper->func, paper->tolerance);
    Roots **roots = paper->roots;

    display_roots(roots, degree);
    return paper;
}


Polynomial *create_equation(unsigned short degree) {
    Polynomial *func = (Polynomial*) malloc(sizeof(Polynomial));
    func->terms = (float *) malloc (degree * sizeof(float));
    func->degree = degree;
    for(int i=0; i<degree;i++) {
        printf("Enter Coefficient (x^%d): ", degree-i-1);
        scanf("%f", &func->terms[i]);
    }


    return func;
}

float delta_u_formula(Polynomial *b_func, Polynomial *c_func) {
    unsigned short max = b_func->degree - 1;

    return ((b_func->terms[max] * c_func->terms[max-3]) - (b_func->terms[max-1] * c_func->terms[max-2])) /
            ((c_func->terms[max-1] * c_func->terms[max-3])-(c_func->terms[max-2] * c_func->terms[max-2]));
}

float delta_v_formula(Polynomial *b_func, Polynomial *c_func) {
    unsigned short max = b_func->degree - 1;

    return  ((c_func->terms[max-1] * b_func->terms[max-1])-(b_func->terms[max] * c_func->terms[max-2])) / 
            ((c_func->terms[max-1] * c_func->terms[max-3])-(c_func->terms[max-2] * c_func->terms[max-2]));

}

Polynomial *equate(float neg_u0, float neg_v0, Polynomial *func) {
    Polynomial *new_func = (Polynomial*)malloc(func->degree * sizeof(Polynomial));
    new_func->degree = func->degree;

    new_func->terms = (float*)malloc(func->degree * sizeof(float));

    for(int i=0;i<func->degree;i++) {
        if(i==0) {
            new_func->terms[i] = func->terms[i];
        } else if(i==1) {
            new_func->terms[i] = func->terms[i] + (neg_u0 * new_func->terms[i-1]);
        }
        else {
            new_func->terms[i] = func->terms[i] + (neg_u0 * new_func->terms[i-1]) + (neg_v0 * new_func->terms[i-2]);
        }
    }

    return new_func;
}

Table *create_table(float u0, float v0, Polynomial *first_func) {

    Table *new_table = (Table*) malloc(sizeof(Table));
    new_table->neg_u0 = u0 * (-1);
    new_table->neg_v0 = v0 * (-1);
    new_table->first_func = first_func;
    new_table->middle_func = equate(new_table->neg_u0, new_table->neg_v0, new_table->first_func);
    new_table->last_func = equate(new_table->neg_u0, new_table->neg_v0, new_table->middle_func);
    display(new_table->first_func);
    display(new_table->middle_func);
    display(new_table->last_func);

    new_table->delta_u = delta_u_formula(new_table->middle_func, new_table->last_func);
    new_table->delta_v = delta_v_formula(new_table->middle_func, new_table->last_func);
    new_table->u_k1 = u0 + new_table->delta_u;
    new_table->v_k1 = v0 + new_table->delta_v;
    printf("U0 = %.6f\tV0=%.6f\n", new_table->delta_u, new_table->delta_v);
    printf("Uk+1 = %.6f\tVk+1 = %.6f\n\n\n", new_table->u_k1, new_table->v_k1);
    return new_table;
}

void display(Polynomial *func) {

    for(int i = 0; i < func->degree; i++) {
        float terms = func->terms[i];
        char op = '+';
        if(func->terms[i] < 0) {
            terms = func->terms[i] * -1;
            op = '-';
        }
        if(i == 0)
            printf("%.6fx^%d", terms, func->degree - i-1);
        else {
            
            printf(" %c %.6fx^%d", op, terms, func->degree - i-1);
        }
            
    }
    printf("\n");

}