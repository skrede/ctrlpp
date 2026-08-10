float promoted(float x)
{
    return static_cast<float>(x * 2.0);
}

int main()
{
    float narrowed = 0.1;
    promoted(narrowed);
    return 0;
}
