# Generated by Django 5.0.3 on 2024-08-09 17:04

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='data',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('refseq_ids', models.CharField(max_length=100)),
                ('variation_ids', models.CharField(max_length=100)),
                ('meanProb', models.FloatField()),
                ('stdProb', models.FloatField()),
                ('pred_label', models.CharField(max_length=100)),
            ],
        ),
    ]
