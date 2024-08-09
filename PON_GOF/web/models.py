from django.db import models


# Create your models here.
class data(models.Model):
    refseq_ids = models.CharField(max_length=100, null=False, blank=False)
    variation_ids = models.CharField(max_length=100, null=False, blank=False)
    meanProb = models.FloatField(null=False, blank=False)
    stdProb = models.FloatField(null=False, blank=False)
    pred_label = models.CharField(max_length=100, null=False, blank=False)

    def __str__(self):
        return self.refseq_ids


class Counter(models.Model):
    val = models.IntegerField(null=False, blank=False, default=0)

    def increment(self):
        self.val += 1
        self.save()

    def __str__(self):
        return str(self.val)
