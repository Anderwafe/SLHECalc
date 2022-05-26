#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

int koefCount;

typedef struct node
{
    float *coefs;
    struct node *next;
} expression;

int addNewExpression(expression **first, float *newCoefs)
{
    expression *newExpression = (expression*)malloc(sizeof(expression));
    if(newExpression == NULL)
    {
        return 1;
    }
    newExpression->coefs = (float*)malloc(sizeof(float)*koefCount);
    for(int i = 0; i < koefCount; i++)
    {
        newExpression->coefs[i] = newCoefs[i];
    }
    newExpression->next = NULL;

    if(*first == NULL)
    {
        *first = newExpression;
        return 0;
    }

    expression *curr = *first;
    while(curr->next != NULL)
    {
        curr = curr->next;
    }

    curr->next = newExpression;


    return 0;
}

expression* getExpression(expression *first, int num)
{
    if(num < 1 || first == NULL) return NULL;

    expression *curr = first;
    for(int i = 1; i < num; i++)
    {
        if(curr->next == NULL)
            return NULL;
        curr = curr->next;
    }
    return curr;
}

expression* expressionSum(expression *addWith, expression *summand)
{
    expression *newExpression = (expression*)malloc(sizeof(expression));
    newExpression->coefs = (float*)malloc(sizeof(float)*koefCount);
    for(int i = 0; i < koefCount; i++)
    {
        newExpression->coefs[i] = addWith->coefs[i] + summand->coefs[i];
    }
    return newExpression;
}

expression* expressionMultiply(expression *expr, float multiplier)
{
    expression *newExpression = (expression*)malloc(sizeof(expression));
    newExpression->coefs = (float*)malloc(sizeof(float)*koefCount);
    for(int i = 0; i < koefCount; i++)
    {
        newExpression->coefs[i] = expr->coefs[i] * multiplier;
    }
    return newExpression;
}

expression* expressionMinor(expression* expr, int row, int column)
{
    row++;
    expression* ret = NULL;
    int count = 1;
    expression* curr = expr;
    while((curr = curr->next) != NULL) count++;
    for(int i = 1; i <= count; i++)
    {
        if(i == row)
        {
            continue;
        }
        float newCoefs[count];
        int j1 = 0;
        for(int j = 0; j < count; j++)
        {
            if(column == j) continue;
            newCoefs[j1++] = getExpression(expr, i)->coefs[j];
        }
        addNewExpression(&ret, newCoefs);
    }
    return ret;
}

float determinantMatrix(expression* expr)
{
    if(expr == NULL) return 100;

    int count = 1;
    expression* curr = expr;
    while((curr = curr->next) != NULL) count++;
    if(count == 1)
    {
        return expr->coefs[0];
    }
    if(count == 2)
    {
        return (((expr->coefs[0])*(expr->next->coefs[1]))-((expr->coefs[1])*(expr->next->coefs[0])));
    }

    float dets[koefCount-1];
    for(int k = 0; k < koefCount-1; k++)
    {
        dets[k] = expr->coefs[k] * determinantMatrix(expressionMinor(expr, 0, k));
    }
    float resDet = 0;
    for(int i = 1; i < koefCount; i++)
    {
        resDet += dets[i-1] * (i % 2 == 0 ? -1 : 1);
    }

    return resDet;
}

float* GauseMethod(expression *expr, int exprCount)
{
    if(expr == NULL)
    {
        printf("expression error 100");
        return 0;
    }
    if(koefCount == 2)
    {
        printf("(x) = %f\n", (expr->coefs[1])/(expr->coefs[0]));

        return expr->coefs;
    }
    for(int i = 1; i < exprCount; i++)
    {
        for(int j = 1; j <= exprCount-i; j++)
        {
            float buf[koefCount];
            for(int k = 0; k < koefCount; k++)
            {
                buf[k] = expressionSum(expressionMultiply(getExpression(expr, i), getExpression(expr, i+j)->coefs[i-1]*(-1)), expressionMultiply(getExpression(expr, i+j), getExpression(expr, i)->coefs[i-1]))->coefs[k];
            }
            for(int k = 0; k < koefCount; k++)
            {
                getExpression(expr, i+j)->coefs[k] = buf[k];
            }
        }
    }

    printf("\n\nС начала приведём систему в треугольный вид.\n");

    expression *curr = expr;
    while(curr != NULL)
    {
        for(int j = 0; j < koefCount-2; j++)
        {
            printf("%f(%d) + ", curr->coefs[j], j+1);
        }
        printf("%f(%d) = %f\n", curr->coefs[koefCount-2], koefCount-1, curr->coefs[koefCount-1]);
        curr = curr->next;
    }

    curr = expr;
    while(curr != NULL)
    {
        for(int j = 0; j < koefCount-1; j++)
        {
            if(curr->coefs[j] != 0)
            {
                curr = curr->next;
                continue;
            }
            if(j == koefCount-1)
            {
                printf("\n\nТак как в системе присутствует уравнение, состоящее из нулей, то система имеет бесконечное количество решений.\n");
                return NULL;
            }
        }
    }

    printf("\n\nПосчитаем, начиная с самого нижнего уравнения системы, все переменные.\n");

    float res[koefCount-1];
    for(int i = exprCount; i > 0; i--)
    {
        float resBuf = getExpression(expr, i)->coefs[koefCount-1];
        for(int j = koefCount-2; j >= i; j--)
        {
            resBuf -= getExpression(expr, i)->coefs[j]*res[j];
        }
        res[i-1] = resBuf/(getExpression(expr, i)->coefs[i-1]);
    }

    for(int j = 0; j < koefCount-1; j++)
    {
        printf("%f(%d); ", res[j], j+1);
    }
    printf("\n");
    return expr->coefs;
}

float* CramerMethod(expression *expr)
{
    printf("\n\nС начала найдем детерминант матрицы, основанной на уравнениях системы.\n");
    float mainDet = determinantMatrix(expr);
    printf("\n\nОпределитель матрицы равен - det(matrix) = %f\n", mainDet);
    printf("\n\nНайдем остальные определители.\n");
    float remainingDet[koefCount-1];
    for(int k = 0; k < koefCount-1; k++)
    {
        expression* remExpr = NULL;
        float remCoefs[koefCount-1];
        for(int i = 1; i < koefCount; i++)
        {
            for(int j = 0; j < koefCount-1; j++)
            {
                remCoefs[j] = (k == j ? getExpression(expr, i)->coefs[koefCount-1] : getExpression(expr, i)->coefs[j]);
            }
            addNewExpression(&remExpr, remCoefs);
        }
        remainingDet[k] = determinantMatrix(remExpr);
        printf("det(%d) = %f\n", k+1, remainingDet[k]);
    }
    printf("\n\nПосчитаем каждую переменную системы по формуле - \n");
    for(int i = 0; i < koefCount-1; i++)
    {
        printf("(%d) = det(%d)/det(matrix) = %f\n", i+1, i+1, remainingDet[i]/mainDet);
    }
    return expr->coefs;
}

float* matrixMethod(expression* expr)
{
    printf("\n\nС начала найдем обратную матрицу.\n");
    int count = 1;
    expression* reverseExpr = NULL;
    expression* curr = expr;
    while((curr = curr->next) != NULL) count++;
    for(int i = 1; i <= count; i++)
    {
        printf("\n| ");
        float reverseParams[count];
        for(int j = 0; j < count; j++)
        {
            reverseParams[j] = determinantMatrix(expressionMinor(expr, j, i-1)) * ((i+j) % 2 == 0 ? -1 : 1) * (1.0/determinantMatrix(expr));
            printf("%15f ", reverseParams[j]);
        }
        addNewExpression(&reverseExpr, reverseParams);
        printf("|");
    }

    printf("\n\nТеперь, перемножим обратную матрицу на свободные элементы СЛАУ. Элементы полученной матрицы - значения неизвестных.\n");

    int countRow = count;

    expression* res = NULL;
    for(int i = 0; i < countRow; i++)
    {
        float newRow[1];
        newRow[0] = 0;
        for(int j = 0; j < countRow; j++)
        {
            newRow[0] += getExpression(reverseExpr, i+1)->coefs[j] * getExpression(expr,j+1)->coefs[koefCount-1];
        }
        printf("(%d) = %f\n", i+1, newRow[0]);
        addNewExpression(&res, newRow);
    }

    return expr->coefs;
}

int main()
{
    expression *head = NULL;
    //float* res;
    start:
    free(head);
    head = NULL;

    //system("clear");
    //res = NULL;
    printf("Введите количество неизвестных системы: ");
    scanf("%d", &koefCount);
    koefCount++;

    system("clear");

    int expressionsCount = 0;
    for(int i = 0; i < koefCount-1; i++)
    {
        if(expressionsCount)
        {
            printf("Добавить еще 1 уравнение?\n1 - да\n2 - нет\n");
            int a;
            scanf("%d", &a);
            if(a == 2) break;
        }
        system("clear");
        float coefs[koefCount];
        for(int i = 0; i < koefCount; i++)
        {
            if(i == koefCount-1)
            {
                printf("Введите свободный член %d уравнения системы: ", expressionsCount+1);
            }
            else
            {
                printf("Введите %d коэффициент %d уравнения системы (учитывайте, что, если переменная в уравнении отсутствует, то коэффициент равен 0): ", i, expressionsCount+1);
            }
            scanf("%f", &(coefs[i]));
        }

        if(addNewExpression(&head, coefs))
        {
            printf("Ошибка выделения памяти, повторите \n");
            i--;

            continue;
        }
        expressionsCount++;
        system("clear");
    }
    if(expressionsCount < koefCount-1)
    {
        printf("Система имеет бесконечное количество решений.\n");
        free(head);
        return 0;
    }

    system("clear");

    while(1)
    {
        system("clear");
        printf("Система имеет вид:\n");
        expression *curr = head;
        while(curr != NULL)
        {
            for(int j = 0; j < koefCount-2; j++)
            {
                printf("%f(%d) + ", curr->coefs[j], j+1);
            }
            printf("%f(%d) = %f\n", curr->coefs[koefCount-2], koefCount-1, curr->coefs[koefCount-1]);
            curr = curr->next;
        }

        printf("Выберите действие:\n1. Решить\n2. Изменить\n3. Начать с начала\n");
        int choice;
        scanf("%d", &choice);
        switch(choice)
        {
            case 1:
            {
                printf("Каким способом её нужно решить?\n1. Метод Гаусса\n2. Метод Крамера\n3. Матричный метод\n");
                int choiceAnother;
                scanf("%d", &choiceAnother);
                switch(choiceAnother)
                {
                    case 1:
                    {
                        GauseMethod(head, expressionsCount);
                        break;
                    }
                    case 2:
                    {
                        CramerMethod(head);
                        break;
                    }
                    case 3:
                    {
                        matrixMethod(head);
                        break;
                    }
                }

                fflush(stdin);
                getchar();

                goto start;
            }
            case 2:
            {
                printf("Какое уравнение вы хотите изменить?\n");
                int changeNum;
                scanf("%d", &changeNum);
                float coefs[koefCount];
                for(int i = 0; i < koefCount; i++)
                {
                    if(i == koefCount-1)
                    {
                        printf("Введите свободный член %d уравнения системы: ", changeNum);
                    }
                    else
                    {
                        printf("Введите %d коэффициент %d уравнения системы (учитывайте, что, если переменная в уравнении отсутствует, то коэффициент равен 0): ", i, changeNum);
                    }
                    scanf("%f", &(coefs[i]));
                }

                for(int k = 0; k < koefCount; k++)
                {
                    getExpression(head, changeNum)->coefs[k] = coefs[k];
                }
                break;
            }
            case 3:
            {
                system("clear");
                goto start;
            }
            default:
            {
                printf("Неверный ввод. Повторите попытку.\n");
                break;
            }
        }
    }
    return 0;
}
